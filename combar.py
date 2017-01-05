#!/usr/bin/python
from combar.combarall import ComBaR
import logging

if __name__ == '__main__':
    try:
        ppw = ComBaR()
        ppw.run()
    except Exception as exception:
        # Show complete exception when running in DEBUG
        if (hasattr(ppw.settings, 'loglevel') and
            getattr(logging, 'DEBUG') == ppw.logger.getEffectiveLevel()):
            ppw.logger.exception(str(exception))
        else:
            print('Program ended. The message was: ', ','.join(exception.args))
            print("Please use the option --help for information on command line arguments.")
