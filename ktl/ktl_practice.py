#! @KPYTHON@

import ktl
import time

from astropy.io import fits
from photometry import photometry
import quadratic


class Keyword:
    def __init__(self, record, pane, exposure, speed):

        self.record = record
        self.pane = pane
        self.exposure = exposure
        self.speed = speed
        self.start = 'Yes'

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

        self.event_key = ktl.cache('nickucam', 'EVENT')
        self.event_key.callback(self.event_callback)
        self.event_key.monitor()
   
    def event_callback(self, keyword):
        self.event_value = keyword.read()
        print(f'update EVENT: {self.event_value}')

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

        print(f'POCSECPD: {self.keyword.secpd_key.read()}')
        print(f'POCSECPA: {self.keyword.secpa_key.read()}')

        if abs(float(self.keyword.secpd_key.read()) - focus_value) < .1:
            print(f'POCSECPD already set to {focus_value}. No change needed.')
            return
        
        if focus_value < 165 or focus_value > 500:
            raise ValueError(f"Focus value {focus_value} is out of range (165-500).")

        if not self.keyword.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot take exposure.")

        print(f'POCSECLK: {self.keyword.seclk_key.read()}')
        self.keyword.seclk_key.write('off')
        print(f'POCSECLK: {self.keyword.seclk_key.read()}')

        self.keyword.secpd_key.write(focus_value)
        print(f'POCSECPD: {self.keyword.secpd_key.read()}')
        print(f'POCSECPA: {self.keyword.secpa_key.read()}')

        if not self.keyword.seclk_key.waitFor('== on', timeout=15):
            raise Exception("POCSECLK did not turn on. Focus change failed.")
    ### CHANGE FOCUS ###

    ### TAKE EXPOSURE ###
    def exposure(self):

        if not self.keyword.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot take exposure.")

        self.keyword.record_key.write(self.keyword.record)
        record_value = self.keyword.record_key.read()
        print(f'RECORD: {record_value}')

        self.keyword.pane_key.write(self.keyword.pane)
        pane_value = self.keyword.pane_key.read()
        print(f'PANE: {pane_value}')

        self.keyword.exposure_key.write(self.keyword.exposure)
        exposure_value = self.keyword.exposure_key.read()
        print(f'EXPOSURE: {exposure_value}')

        self.keyword.speed_key.write(self.keyword.speed)
        speed_value = self.keyword.speed_key.read()
        print(f'READSPEED: {speed_value}')

        self.keyword.start_key.write(self.keyword.start)
        start_value = self.keyword.start_key.read()
        print(f'START: {start_value}')

        self.keyword.event_key.waitFor('== ReadoutBegin', timeout=round(self.keyword.exposure * 1.2))

        if self.keyword.event_key.waitFor('== ReadoutEnd', timeout=30):
            pass
    ### TAKE EXPOSURE ###

    def sequence(self, focus_value):

        filepath = f"{self.keyword.dir_key.read()}/{self.keyword.file_key.read()}{self.keyword.obs_key.read()}.{self.keyword.suffix_key.read()}"
        print(f"Exposure being saved at: {filepath}")

        # if self.keyword.record == 'Yes':
        #     self.keyword.dir_key.write("/data")
        #     self.keyword.file_key.write(f"{count}focus_")
        #     self.keyword.obs_key.write(str(focus_value))
        #     self.keyword.suffix_key.write("fits")

        self.focus(focus_value)
        self.focus_value = focus_value
        self.exposure()

        # hdu = fits.open(filepath)
        # print(hdu.info())

        # self.fwhm = photometry(filepath, verbose=False)
        # print(f" {self.keyword.file_key.read()}{self.keyword.obs_key.read()}.{self.keyword.suffix_key.read()} FWHM: {fwhm} \n")

def curve_finder(image1, image2, seen, keyword, direction=None):
    seen.add(image1)
    seen.add(image2)
    
    if image1.fwhm > image2.fwhm:
        if direction is None:
            direction = 'right'
        if direction == 'left':
            return curve_helper(image1, image2, seen)
        focus3 = image2.focus_value + (image2.focus_value - image1.focus_value)
        image3 = Event(keyword)
        image3.sequence(focus3)
        return curve_finder(image2, image3, seen, direction)
    elif image1.fwhm < image2.fwhm:
        if direction is None:
            direction = 'left'
        if direction == 'right':
            return curve_helper(image1, image2, seen)
        focus3 = image1.focus_value - (image2.focus_value - image1.focus_value)
        image3 = Event(keyword)
        image3.sequence(focus3)
        return curve_finder(image3, image1, seen, direction)
    
def curve_helper(image1, image2, seen, iterations=2, keyword):
    if iterations > 0:
        if image1.fwhm > image2.fwhm:
            focus3 = image1.focus_value - (image2.focus_value - image1.focus_value)
            image3 = Event(keyword)
            image3.sequence(focus3)
            seen.add(image3)
            return curve_helper(image3, image1, seen, iterations-1)
        elif image1.fwhm < image2.fwhm:
            focus3 = image2.focus_value + (image2.focus_value - image1.focus_value)
            image3 = Event(keyword)
            image3.sequence(focus3)
            seen.add(image3)
            return curve_helper(image2, image3, seen, iterations-1)
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

    filepath = f"{keyword.dir_key.read()}/{keyword.file_key.read()}{keyword.obs_key.read()}.{keyword.suffix_key.read()}"

    if not keyword.wait_until(keyword.event_key, ['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")


    print("taking exposure")
    event.exposure()
    print("Exposure taken")

    # event.focus(365)

    # event.sequence(360)
    # event.sequence(365)

    # event.sequence(args.focus_value)

    # focus_finder(args.focus_value, 5, keyword)


print("hello world")
# if __name__ == "__main__":
#     main()


# /data/nickel
# 7:30

secpa_key = ktl.cache('nickelpoco', 'POCSECPA')
print(f'secpa: {secpa_key.read()}')