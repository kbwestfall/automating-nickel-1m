#! @KPYTHON@

import ktl
import time

class Keyword:
    def __init__(self):
        self.general = None

        self.stop_key = ktl.cache("nickelpoco", "POCSTOP")

        self.event_key = ktl.cache("nickucam", "EVENT")
        self.event_key.callback(self.event_callback)
        # self.event_key.monitor()

        self.seclk_key = ktl.cache("nickelpoco", "POCSECLK")
        self.secpd_key = ktl.cache("nickelpoco", "POCSECPD")

        self.seclk_key.callback(self.seclk_callback)
        self.seclk_key.monitor()

        self.record_key = ktl.cache('nickucam', 'RECORD')
        self.exposure_key = ktl.cache('nickucam', 'EXPOSURE')
        self.speed_key = ktl.cache('nickucam', 'READSPEED')
        self.start_key = ktl.cache('nickucam', 'START')

        self.record_key.write('off')
        self.exposure_key.write(1)
        self.speed_key.write('fast')

    def general_callback(self, keyword):
        self.general = keyword.read()

    def event_callback(self, keyword):
        self.event = keyword.read()

    def seclk_callback(self, keyword):
        self.seclk = keyword.read()

    def wait_for(self, keyword, expected_value, timeout=15):
        keyword.callback(self.general_callback)
        keyword.monitor(start=True, prime=True)
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self.general in expected_value:
                return True
        return False

import argparse

def main():
    parser = argparse.ArgumentParser(description='waitFor behavior test.')
    parser.add_argument('-f', '--focus', type=int, default=360, help='desired focus value.')
    args = parser.parse_args()

    keyword = Keyword()

    if not keyword.wait_for(keyword.stop_key, ['allowed'], timeout=15):
        error_msg = "POCSTOP not allowed. Cannot change focus."
        raise Exception(error_msg)
    
    if not keyword.wait_for(keyword.event_key, ['ControllerReady', 'ExposeSequenceDone'], timeout=15):
        error_msg = "Controller not ready. Cannot change focus."
        raise Exception(error_msg)

    keyword.seclk_key.write('off')
    keyword.seclk_key.read()
    keyword.secpd_key.write(args.focus)

    # # wait_for
    # if not keyword.wait_for(keyword.seclk_key, ["on"], timeout=60):
    #     error_msg = "SECLK did not turn on."
    #     raise Exception(error_msg)
    # else:
    #     keyword.start_key.write('on')

    #     if not keyword.wait_for(keyword.event_key, ['ExposureBegin'], timeout=60):
    #         error_msg = "Exposure sequence did not start."
    #         raise Exception(error_msg)

    #     if not keyword.wait_for(keyword.event_key, ['ControllerReady'], timeout=60):
    #         error_msg = "Exposure sequence did not complete."
    #         raise Exception(error_msg)
    #     else:
    #         print("Exposure completed")

    # waitFor
    
    if not keyword.seclk_key.waitFor('== on', timeout=60):
        error_msg = "SECLK did not turn on."
        raise Exception(error_msg)
    else:
        keyword.start_key.write('on')

        if not keyword.event_key.waitFor('== ExposureBegin', timeout=60):
            error_msg = "Exposure sequence did not start."
            raise Exception(error_msg)

        if not keyword.event_key.waitFor('== ControllerReady', timeout=60):
            error_msg = "Exposure sequence did not complete."
            raise Exception(error_msg)
        else:
            print("Exposure completed")

if __name__ == "__main__":
    main()