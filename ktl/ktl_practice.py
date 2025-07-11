#! @KPYTHON@

import ktl
import time

from photometry import photometry
from astropy.io import fits


class Event:
    def __init__(self):

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
        self.event_key.callback(self.callback)
        self.event_key.monitor()
        self.event_value = self.event_key.read()

    def callback(self, keyword):
        self.event_value = keyword.read()
        print(f'update EVENT: {self.event_value}')

    def wait_until(self, keyword, expected_value, timeout=15):
        start_time = time.time()
        while time.time() - start_time < timeout:
            if keyword.read() in expected_value:
                return True
            time.sleep(1)
        return False

    ### CHANGE FOCUS ###
    def focus(focus_value):

        print(f'POCSECPD: {self.secpd_key.read()}')
        print(f'POCSECPA: {self.secpa_key.read()}')

        if abs(float(self.secpd_key.read()) - focus_value) < .1:
            print(f'POCSECPD already set to {focus_value}. No change needed.')
            return
        
        if focus_value < 165 or focus_value > 500:
            raise ValueError(f"Focus value {focus_value} is out of range (165-500).")

        if not self.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot take exposure.")

        print(f'POCSECLK: {self.seclk_key.read()}')
        self.seclk_key.write('off')
        print(f'POCSECLK: {self.seclk_key.read()}')

        self.secpd_key.write(focus_value)
        print(f'POCSECPD: {self.secpd_key.read()}')
        print(f'POCSECPA: {self.secpa_key.read()}')

        if not self.seclk_key.waitFor('== on', timeout=15):
            raise Exception("POCSECLK did not turn on. Focus change failed.")
    ### CHANGE FOCUS ###

    ### TAKE EXPOSURE ###
    def exposure(record, window, exposure, speed, start):

        if not self.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot take exposure.")

        self.record_key.write(record)
        record_value = self.record_key.read()
        print(f'RECORD: {record_value}')

        self.pane_key.write(window)
        pane_value = self.pane_key.read()
        print(f'PANE: {pane_value}')

        self.exposure_key.write(exposure)
        exposure_value = self.exposure_key.read()
        print(f'EXPOSURE: {exposure_value}')

        self.speed_key.write(speed)
        speed_value = self.speed_key.read()
        print(f'READSPEED: {speed_value}')

        # set the start keyword to 'Yes'
        self.start_key.write(start)
        start_value = self.start_key.read()
        print(f'START: {start_value}')

        self.event_key.waitFor('== ReadoutBegin', timeout=round(exposure * 1.2))

        if self.event_key.waitFor('== ReadoutEnd', timeout=30):
            pass
    ### TAKE EXPOSURE ###

    def sequence(focus_value, record, window, exposure_length, speed, start):

        filepath = f"{self.dir_key.read()}/{self.file_key.read()}{self.obs_key.read()}.{self.suffix_key.read()}"
        print(f"Exposure being saved at: {filepath}")

        self.focus(focus_value)
        self.exposure(record, window, exposure_length, speed, start)

        # hdu = fits.open(filepath)
        # print(hdu.info())

        # fwhm = photometry(filepath, verbose=False)
        # print(f" {file_key.read()}{obs_key.read()}.{suffix_key.read()} FWHM: {fwhm} \n")

import argparse

def main():

    parser = argparse.ArgumentParser(description="Automate the focus finding process.")
    parser.add_argument('-f', '--focus_value', default=360, help='Focus value')
    parser.add_argument('-e', '--exposure_length', default=1, help='Exposure length in seconds')
    parser.add_argument('-w', '--window_size', default="0 0 2048 2048", help='Window size for exposure')
    args = parser.parse_args()

    event = Event()
    print(f'EVENT: {event.event_value}')

    filepath = f"{event.dir_key.read()}/{event.file_key.read()}{event.obs_key.read()}.{event.suffix_key.read()}"

    if not event.wait_until(['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")


    # print("taking exposure")
    # event.exposure('No', 0, 'Slow', 1)
    # print("Exposure taken")

    # event.focus(365)

    event.sequence(args.focus_value, 'No', args.window_size, args.exposure_length, 'Fast', 'Yes')


    # focus_value = 350
    # while focus_value < 375:
    #     event.sequence(focus_value, 'No', 1, 'Fast', 'Yes')
    #     focus_value += 5


print("hello world")
# if __name__ == "__main__":
#     main()

