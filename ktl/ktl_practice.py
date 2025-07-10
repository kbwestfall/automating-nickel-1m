#! @KPYTHON@

import ktl
import time

# service = ktl.Service('nickelpoco')
# service = ktl.cache('nickelpoco')

# keyword = ktl.cache('nickelpoco', 'POCRAA')
# value = keyword.read()
# print(f'POCRAA: {value}')

# def callback(keyword):
#     new_value = keyword.read()
#     print(f'update POCRAA: {new_value}')

# keyword.callback(callback)
# keyword.monitor()
# time.sleep(1)
# keyword.monitor(start=False)

# secpa = ktl.cache('nickelpoco', 'POCSECPA')
# secpd = ktl.cache('nickelpoco', 'POCSECPD')
# print(f'POCSECPA: {secpa.read()}')
# print(f'POCSECPD: {secpd.read()}')


### TAKE EXPOSURE ###
def exposure(record, window, exposure, speed, start, event):
    record_key = ktl.cache('nickucam', 'RECORD')
    exposure_key = ktl.cache('nickucam', 'EXPOSURE')
    speed_key = ktl.cache('nickucam', 'READSPEED')
    start_key = ktl.cache('nickucam', 'START')
    pane_key = ktl.cache('nickucam', 'PANE')

    if not event.event_key.waitFor('== ControllerReady', timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")

    record_key.write(record) 
    record_value = record_key.read()
    print(f'RECORD: {record_value}')

    pane_key.write(window)
    pane_value = pane_key.read()
    print(f'PANE: {pane_value}')

    exposure_key.write(exposure)
    exposure_value = exposure_key.read()
    print(f'EXPOSURE: {exposure_value}')

    speed_key.write(speed)
    speed_value = speed_key.read()
    print(f'READSPEED: {speed_value}')

    # set the start keyword to 'Yes'
    start_key.write(start)
    start_value = start_key.read()
    print(f'START: {start_value}')

    event.event_key.waitFor('== ReadoutBegin', timeout=round(exposure * 1.2))

    if event.event_key.waitFor('== ReadoutEnd', timeout=30):
        pass
### TAKE EXPOSURE ###


### CHANGE FOCUS ###
def focus(focus_value, event):
    secpa_key = ktl.cache('nickelpoco', 'POCSECPA')
    secpd_key = ktl.cache('nickelpoco', 'POCSECPD')
    seclk_key = ktl.cache('nickelpoco', 'POCSECLK')

    print(f'POCSECPD: {secpd_key.read()}')
    print(f'POCSECPA: {secpa_key.read()}')

    if abs(float(secpd_key.read()) - focus_value) < .1:
        print(f'POCSECPD already set to {focus_value}. No change needed.')
        return
    
    if focus_value < 165 or focus_value > 500:
        raise ValueError(f"Focus value {focus_value} is out of range (165-500).")

    if not event.event_key.waitFor('== ControllerReady', timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")

    print(f'POCSECLK: {seclk_key.read()}')
    seclk_key.write('off')
    print(f'POCSECLK: {seclk_key.read()}')

    print(f'POCSECPD: {secpd_key.read()}')
    print(f'POCSECPA: {secpa_key.read()}')
    secpd_key.write(focus_value)
    print(f'POCSECPD: {secpd_key.read()}')
    print(f'POCSECPA: {secpa_key.read()}')

    if not seclk_key.waitFor('== on', timeout=15):
        raise Exception("POCSECLK did not turn on. Focus change failed.")
### CHANGE FOCUS ###


### race condition ###
from photometry import photometry
from astropy.io import fits

def sequence(focus_value, record, window, exposure_length, speed, start, event):

    obs_key = ktl.cache('nickucam', 'OBSNUM')
    dir_key = ktl.cache('nickucam', 'OUTDIR')
    file_key = ktl.cache('nickucam', 'OUTFILE')
    suffix_key = ktl.cache('nickucam', 'RECORD_SUFFIX')

    filepath = f"{dir_key.read()}/{file_key.read()}{obs_key.read()}.{suffix_key.read()}"
    print(f"Exposure being saved at: {filepath}")

    focus(focus_value, event)
    exposure(record, window, exposure_length, speed, start, event)

    # hdu = fits.open(filepath)
    # print(hdu.info())

    # fwhm = photometry(filepath, verbose=False)
    # print(f" {file_key.read()}{obs_key.read()}.{suffix_key.read()} FWHM: {fwhm} \n")
### race condition ###


class Event:
    def __init__(self):
        self.event_key = ktl.cache('nickucam', 'EVENT')
        self.event_key.callback(self.callback)
        self.event_key.monitor()
        self.event_value = self.event_key.read()

    def callback(self, keyword):
        self.event_value = keyword.read()
        print(f'update EVENT: {self.event_value}')

    def wait_for(self, expected_value, timeout=15):
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self.event_value in expected_value:
                return True
            time.sleep(1)
        return False

import argparse

def main():

    parser = argparse.ArgumentParser(description="Automate the focus finding process.")
    parser.add_argument('-f', '--focus_value', default=360, help='Focus value')
    parser.add_argument('-e', '--exposure_length', default=1, help='Exposure length in seconds')
    parser.add_argument('-w', '--window_size', default="0 0 2048 2048", help='Window size for exposure')
    args = parser.parse_args()

    event = Event()
    print(f'EVENT: {event.event_value}')

    if not event.wait_for(['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")


    # print("taking exposure")
    # exposure('No', 0, 'Slow', 1, event)
    # print("Exposure taken")

    #focus(365, event)

    sequence(args.focus_value, 'No', args.window_size, args.exposure_length, 'Fast', 'Yes', event)


    # focus_value = 350
    # while focus_value < 375:
    #     sequence(focus_value, 'No', 1, 'Fast', 'Yes', event)
    #     focus_value += 5


print("hello world")
# if __name__ == "__main__":
#     main()

