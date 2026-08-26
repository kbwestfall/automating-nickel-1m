#! @KPYTHON@

import ktl
import time

class Keyword:
    def __init__(self):
        self.general = None

        self.integer_key = ktl.cache("pie", "INTEGER")
        self.integer_key.callback(self.integer_callback)
        self.integer_key.monitor()

        self.boolean_key = ktl.cache("pie", "BOOLEAN")
        self.boolean_key.callback(self.boolean_callback)
        self.boolean_key.monitor()

    def general_callback(self, keyword):
        print(f'GENERAL CHANGE: {keyword.name} = {keyword.read()}')
        self.general = keyword.read()

    def boolean_callback(self, keyword):
        print(f'BOOLEAN CHANGE: {keyword.read()}')

    def integer_callback(self, keyword):
        print(f'INTEGER CHANGE: {keyword.read()}')

    def wait_for(self, keyword, expected_value, timeout=15):
        keyword.callback(self.general_callback)
        keyword.monitor(start=True)
        keyword.propagate()
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self.general in expected_value:
                return True
        return False

def main():
    keyword = Keyword()
    keyword.integer_key.write(42)
    
    if not keyword.wait_for(keyword.integer_key, [100], timeout=60):
        keyword.integer_key.write(500)

    if not keyword.wait_for(keyword.boolean_key, ["True"], timeout=60):
        keyword.integer_key.write(0)

    # if not keyword.boolean_key.waitFor('== True', timeout=50):
    #     keyword.integer_key.write(0)

    else:
        keyword.integer_key.write(1)

if __name__ == "__main__":
    main()