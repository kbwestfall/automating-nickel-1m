#! @KPYTHON@

import ktl
import time

from astropy.io import fits
from photometry import photometry
from simulation import simulation
import quadratic

class Image:
    def __init__(self, focus):
        self.focus = focus

        true_fwhm = 0.016*(focus-360)**2+6.0

        file_name = str(true_fwhm)
        self.filepath = f"data/nickel/focus{file_name.replace('.', '_')}"

        simulation(self.filepath, fwhm=true_fwhm)


class Keyword:
    def __init__(self, record, pane, exposure, speed):

        self.record = record
        self.pane = pane
        self.exposure = exposure
        self.speed = speed
        self.start = 'No'

        self.secpa_key = ktl.cache('nickelpoco', 'POCSECPA')
        self.secpd_key = ktl.cache('nickelpoco', 'POCSECPD')
        self.seclk_key = ktl.cache('nickelpoco', 'POCSECLK')
        self.obs_key = ktl.cache('nickucam', 'OBSNUM')
        self.dir_key = ktl.cache('nickucam', 'OUTDIR')
        self.file_key = ktl.cache('nickucam', 'OUTFILE')
        self.suffix_key = ktl.cache('nickucam', 'RECORD_SUFFIX')
        self.record_key = ktl.cache('nickucam', 'RECORD')
        self.exposure_key = ktl.cache('nickucam', 'EXPOSURE')
        self.speed_key = ktl.cache('nickucam', 'READSPEED')
        self.start_key = ktl.cache('nickucam', 'START')
        self.pane_key = ktl.cache('nickucam', 'PANE')

        self.stop_key = ktl.cache('nickelpoco', 'POCSTOP')

        self.event_key = ktl.cache('nickucam', 'EVENT')
        self.event_key.callback(self.event_callback)
        self.event_key.monitor()

        self.start_key.callback(self.filepath_callback)
        self.start_key.monitor()

    def event_callback(self, keyword):
        self.event_value = keyword.read()
        print(f'update EVENT: {self.event_value}')

    def filepath_callback(self, keyword):
        self.filepath = f"{self.dir_key.read()}/nickel/{self.file_key.read()}{self.obs_key.read()}.{self.suffix_key.read()}"
        print(f'update FILEPATH: {self.filepath}')

    def wait_until(self, keyword, expected_value, timeout=15):
        start_time = time.time()
        while time.time() - start_time < timeout:
            if keyword.read() in expected_value:
                return True
            time.sleep(1)
        return False

        

class Event:
    def __init__(self, keyword):

        self.keyword = keyword

    ### CHANGE FOCUS ###
    def focus(self, focus_value):

        if abs(float(self.keyword.secpd_key.read()) - focus_value) < .1:
            print(f'POCSECPD already set to {focus_value}. No change needed.')
            return
        
        if focus_value < 165 or focus_value > 500:
            raise ValueError(f"Focus value {focus_value} is out of range (165-500).")

        if not self.keyword.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot take exposure.")

        print(f'POCSECLK: {self.keyword.seclk_key.read()}')
        # set POCSECLK to 'off'

        # change POCSECPD to desired focus value
        print(f'POCSECPD: {self.keyword.secpd_key.read()}')
        print(f'POCSECPA: {self.keyword.secpa_key.read()}')

        if not self.keyword.seclk_key.waitFor('== on', timeout=15):
            raise Exception("POCSECLK did not turn on. Focus change failed.")
    ### CHANGE FOCUS ###

    ### TAKE EXPOSURE ###
    def exposure(self):

        if not self.keyword.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot take exposure.")

        
        record_value = self.keyword.record_key.read()
        print(f'RECORD: {record_value}')

        
        pane_value = self.keyword.pane_key.read()
        print(f'PANE: {pane_value}')

        
        exposure_value = self.keyword.exposure_key.read()
        print(f'EXPOSURE: {exposure_value}')

        
        speed_value = self.keyword.speed_key.read()
        print(f'READSPEED: {speed_value}')

        
        start_value = self.keyword.start_key.read()
        print(f'START: {start_value}')

        self.keyword.event_key.waitFor('== ReadoutBegin', timeout=round(self.keyword.exposure * 1.2))

        if self.keyword.event_key.waitFor('== ReadoutEnd', timeout=30):
            pass
    ### TAKE EXPOSURE ###

    def sequence(self, focus_value):

        # if self.keyword.record == 'Yes':
        #     self.keyword.dir_key.write("/data")
        #     self.keyword.file_key.write(f"{count}focus_")
        #     self.keyword.obs_key.write(str(focus_value))
        #     self.keyword.suffix_key.write("fits")

        # not needed for simulation
        # self.focus(focus_value)
        self.focus_value = focus_value
        # self.exposure()

        print(f'Focus value: {focus_value}')

        img = Image(focus_value)

        filepath = img.filepath
        # print(f"Exposure being saved at: {filepath}")

        self.fwhm = photometry(f"{filepath}.fits", verbose=True)
        print(f" {filepath} FWHM: {self.fwhm} \n")


def curve_finder(image1, image2, seen, keyword, direction=None):
    seen.add(image1)
    seen.add(image2)
    print(f"Curve finder: image1: {image1.focus_value}, FWHM: {image1.fwhm} | image2: {image2.focus_value}, FWHM: {image2.fwhm} | direction: {direction}\n")

    if image1.fwhm > image2.fwhm:
        if direction is None:
            direction = 'right'
        if direction == 'left':
            return curve_helper(image1, image2, seen, keyword)
        focus3 = image2.focus_value + (image2.focus_value - image1.focus_value)
        image3 = Event(keyword)
        image3.sequence(focus3)
        return curve_finder(image2, image3, seen, keyword, direction)
    elif image1.fwhm < image2.fwhm:
        if direction is None:
            direction = 'left'
        if direction == 'right':
            return curve_helper(image1, image2, seen, keyword)
        focus3 = image1.focus_value - (image2.focus_value - image1.focus_value)
        image3 = Event(keyword)
        image3.sequence(focus3)
        return curve_finder(image3, image1, seen, keyword, direction)

def curve_helper(image1, image2, seen, keyword, iterations=2):
    print(f"Curve helper: image1: {image1.focus_value}, FWHM: {image1.fwhm} | image2: {image2.focus_value}, FWHM: {image2.fwhm} \n")

    if iterations > 0:
        if image1.fwhm > image2.fwhm:
            focus3 = image1.focus_value - (image2.focus_value - image1.focus_value)
            image3 = Event(keyword)
            image3.sequence(focus3)
            seen.add(image3)
            return curve_helper(image3, image1, seen, keyword,iterations-1)
        elif image1.fwhm < image2.fwhm:
            focus3 = image2.focus_value + (image2.focus_value - image1.focus_value)
            image3 = Event(keyword)
            image3.sequence(focus3)
            seen.add(image3)
            return curve_helper(image2, image3, seen, keyword,iterations-1)
    else:
        return seen

def focus_finder(initial_focus, step_size, keyword):

    seen = set()
    
    image1 = Event(keyword)
    image1.sequence(initial_focus)

    image2 = Event(keyword)
    image2.sequence(initial_focus + step_size)

    curve = curve_finder(image1, image2, seen, keyword)
    curve = sorted(curve, key=lambda img: img.focus_value)
    x_values = []
    y_values = []
    print("Curve found with the following focus values:")
    for img in curve:
        print(f"Focus: {img.focus_value}, FWHM: {img.fwhm}")
        x_values.append(img.focus_value)
        y_values.append(img.fwhm)

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")


import argparse

def main():

    parser = argparse.ArgumentParser(description="Automate the focus finding process.")
    parser.add_argument('-f', '--focus_value', default=360, help='Focus value')
    parser.add_argument('-e', '--exposure_length', default=1, help='Exposure length in seconds')
    parser.add_argument('-w', '--window_size', default="0 0 2048 2048", help='Window size for exposure')
    args = parser.parse_args()

    keyword = Keyword('No', args.window_size, args.exposure_length, 'Fast')

    event = Event(keyword)
    print(f'EVENT: {keyword.event_key.read()}')

    filepath = f"{keyword.dir_key.read()}/nickel/{keyword.file_key.read()}{keyword.obs_key.read()}.{keyword.suffix_key.read()}"

    if not keyword.wait_until(keyword.event_key, ['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")

    #check if pocstop is enabled
    if keyword.stop_key.read() != 'enabled':
        print("POCSTOP is 'disabled'. Waiting for 'enabled' to allow motion")
    # POCSTOP not needed for simulation
    # if not keyword.stop_key.waitFor('== 0', timeout=10):
    #     raise Exception("POCSTOP is 'disabled'. Set to 'enabled' to allow motion")

    focus_finder(360, 5, keyword)



if __name__ == "__main__":
    main()


# /data/nickel
# 7:30


# find good keyword to monitor for filename callback

# work on visual interface aspect; show focus star with centroid and pane, focus curve