#! @KPYTHON@

import ktl

import time

from photometry import photometry
from focus import Image
import quadratic

class Exposure:
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

        self.event_key = ktl.cache('nickucam', 'EVENT')
        self.event_key.callback(self.event_callback)
        self.event_key.monitor()
        self.event_value = self.event_key.read()
    
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

    def change_focus(self, focus_value):
        print(f'POCSECPA: {self.secpa_key.read()}')
        print(f'POCSECPD: {self.secpd_key.read()}')
        
        if abs(float(self.secpd_key.read()) - float(focus_value)) < 0.1:
            print(f"Focus already set to {focus_value}. No change needed.")
            return

        if focus_value < 165 or focus_value > 500:
            raise ValueError(f"Focus value {focus_value} is out of bounds (165-500).")
        
        if not self.event_key.waitFor('== ControllerReady', timeout=15):
            raise Exception("Controller not ready. Cannot change focus.")

        # set POCSECLOK to 'off'
        # set POCSECPD to value specified by user
        # set POCSECLOK to 'on'
        # print(f'POCSECLK: {self.seclk_key.read()}')
        self.focus = self.secpd_key.read()
        print(f'Changing focus to: {focus_value}')

        ##### FOR TESTING PURPOSES ONLY #####
        self.focus = focus_value
        ##### FOR TESTING PURPOSES ONLY #####

    def take_exposure(self):

        if not wait_until(self.event_key, 'ControllerReady'):
            raise Exception("Controller not ready. Cannot take exposure.")

        obs_num = self.obs_key.read()
        # print(f'OBS: {obs_num}')
        out_dir = self.dir_key.read()
        # print(f'DIR: {out_dir}')
        out_file = self.file_key.read()
        # print(f'FILE: {out_file}')
        suffix = self.suffix_key.read()
        # print(f'SUFFIX: {suffix}')
        filepath = f"{out_dir}/{out_file}{obs_num}.{suffix}"
        print(f"Exposure being saved at: {filepath}")

        # set the record keyword to 'Yes'
        record_value = self.record_key.read()
        print(f'RECORD: {record_value}')

        exposure_value = self.exposure_key.read()
        print(f'EXPOSURE: {exposure_value}')

        # set the read speed to 'fast'
        speed_value = self.speed_key.read()
        print(f'READSPEED: {speed_value}')

        # set the start keyword to 'Yes'
        start_value = self.start_key.read()
        print(f'START: {start_value}')

        ##### FOR TESTING PURPOSES ONLY #####
        Image(self.focus)
        true_fwhm = 0.016*(self.focus-360)**2+6.0
        file_name = str(true_fwhm)
        file_name = file_name.replace('.', '_')
        filepath = f'images/focus_{file_name}.fits'
        ##### FOR TESTING PURPOSES ONLY #####

        self.fwhm = photometry(filepath, verbose=False)
        print(f" {out_file}{obs_num}.{suffix} FWHM: {self.fwhm} \n")


def curve_finder(image1, image2, seen, direction=None):
    seen.add(image1)
    seen.add(image2)
    
    if image1.fwhm > image2.fwhm:
        if direction is None:
            direction = 'right'
        if direction == 'left':
            return curve_helper(image1, image2, seen)
        focus3 = image2.focus + (image2.focus - image1.focus)
        image3 = Exposure()
        image3.change_focus(focus3)
        image3.take_exposure()
        seen.add(image3)
        return curve_finder(image2, image3, seen, direction)
    elif image1.fwhm < image2.fwhm:
        if direction is None:
            direction = 'left'
        if direction == 'right':
            return curve_helper(image1, image2, seen)
        focus3 = image1.focus - (image2.focus - image1.focus)
        image3 = Exposure()
        image3.change_focus(focus3)
        image3.take_exposure()
        seen.add(image3)
        return curve_finder(image3, image1, seen, direction)
    
def curve_helper(image1, image2, seen, iterations=2):
    if iterations > 0:
        if image1.fwhm > image2.fwhm:
            focus3 = image1.focus - (image2.focus - image1.focus)
            image3 = Exposure()
            image3.change_focus(focus3)
            image3.take_exposure()
            seen.add(image3)
            return curve_helper(image3, image1, seen, iterations-1)
        elif image1.fwhm < image2.fwhm:
            focus3 = image2.focus + (image2.focus - image1.focus)
            image3 = Exposure()
            image3.change_focus(focus3)
            image3.take_exposure()
            seen.add(image3)
            return curve_helper(image2, image3, seen, iterations-1)
    else:
        return seen

def take_sequence(initial_focus, step_size=5):

    exp1 = Exposure()
    exp1.change_focus(initial_focus)
    exp1.take_exposure()

    exp2 = Exposure()
    exp2.change_focus(initial_focus + step_size)
    exp2.take_exposure()

    curve = curve_finder(exp1, exp2, seen=set())
    curve = sorted(curve, key=lambda exp: exp.focus)
    x_values = []
    y_values = []
    print("Curve found with the following focus values:")
    for exp in curve:
        print(f"Focus: {exp.focus}, FWHM: {exp.fwhm}")
        x_values.append(exp.focus)
        y_values.append(exp.fwhm)

    a, b, c = quadratic.fit_quadratic(x_values, y_values)
    x_vertex, y_vertex = quadratic.vertex(a, b, c)
    print(f"\nOptimal focus: {x_vertex}, FWHM: {y_vertex}")


def event_callback(keyword):
    event = keyword.read()
    if event != 'ControllerReady':
        # If the event is not 'ControllerReady', we raise an exception
        raise Exception(f"Unexpected event: {event}")

    # If we reach this point, the event is 'ControllerReady'
    # print("Ready to take exposure")
    return

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

def automate(initial_focus, step_size):
    poco = ktl.Service('nickelpoco')
    poco = ktl.cache('nickelpoco')
    ucam = ktl.Service('nickucam')
    ucam = ktl.cache('nickucam')

    event_key = ktl.cache('nickucam', 'EVENT')
    print(f'EVENT: {event_key.read()}')

    event = Event()
    if not event.wait_for(['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        raise Exception("Controller not ready. Cannot take exposure.")

    event_key.callback(event_callback)
    event_key.monitor()

    take_sequence(initial_focus=initial_focus, step_size=step_size)

    

import argparse

def main():
    parser = argparse.ArgumentParser(description="Automate the focus finding process.")
    parser.add_argument('--initial_focus', type=int, default=360, help='Initial focus value to start sequence at')
        # maybe add something that says: if intial focus is larger that 360 then subtract step size, vice versa. 
    parser.add_argument('--step_size', type=int, default=5, help='Step size for focus increments')
    args = parser.parse_args()
    
    automate(initial_focus=args.initial_focus, step_size=args.step_size)


if __name__ == "__main__":
    main()