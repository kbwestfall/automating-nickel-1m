#! @KPYTHON@

import ktl
import time

from photometry import photometry
from astropy.io import fits

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
    def focus(focus_value):

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
    def exposure():

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

    def sequence(focus_value):

        filepath = f"{self.keyword.dir_key.read()}/{self.keyword.file_key.read()}{self.keyword.obs_key.read()}.{self.keyword.suffix_key.read()}"
        print(f"Exposure being saved at: {filepath}")

        self.focus(focus_value)
        self.exposure()

        # hdu = fits.open(filepath)
        # print(hdu.info())

        # fwhm = photometry(filepath, verbose=False)
        # print(f" {file_key.read()}{obs_key.read()}.{suffix_key.read()} FWHM: {fwhm} \n")

def pseudo_focus_finder(initial_focus, step_size, keyword):
    
    images = []
    count = 0
    focus_value = initial_focus
    while focus_value < 375:
        if keyword.record == 'Yes':
            keyword.dir_key.write("/data")
            keyword.file_key.write(f"{count}focus_")
            keyword.obs_key.write(str(focus_value))
            keyword.suffix_key.write("fits")
        images[count] = Event(keyword)
        images[count].sequence(focus_value)
        focus_value += step_size
        count += 1
    return images



import argparse

def main():

    parser = argparse.ArgumentParser(description="Automate the focus finding process.")
    parser.add_argument('-f', '--focus_value', default=360, help='Focus value')
    parser.add_argument('-e', '--exposure_length', default=1, help='Exposure length in seconds')
    parser.add_argument('-w', '--window_size', default="0 0 2048 2048", help='Window size for exposure')
    args = parser.parse_args()

    keyword = Keyword('No', args.window_size, args.exposure_length, 'Fast')

    event = Event(keyword)
    print(f'EVENT: {keyword.event_value}')

    filepath = f"{keyword.dir_key.read()}/{keyword.file_key.read()}{keyword.obs_key.read()}.{keyword.suffix_key.read()}"

    if not keyword.wait_until(['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")


    print("taking exposure")
    event.exposure()
    print("Exposure taken")

    # event.focus(365)

    # event.sequence(360)
    # event.sequence(365)

    # event.sequence(args.focus_value)

    # images = pseudo_focus_finder(float(args.focus_value), 5)
    # print(f"Images taken: {len(images)}")


print("hello world")
# if __name__ == "__main__":
#     main()

