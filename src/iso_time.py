
import numpy as np
import matplotlib.pyplot as plt


class iso_time(object):

    def __init__(self,
                 final_time,
                 delta_time,
                 time_units='seconds',  # 'seconds' , 'minutes', 'hours', 'days'
                 start_time=0):

        self.__time_unit = time_units
        self.__start_time = start_time
        self.__final_time = final_time
        self.dt = delta_time
        self.time_steps = int((self.__final_time - self.__start_time) / self.dt)

    def total_seconds(self):

        seconds = {'seconds': 1, 'minutes': 60, 'hours': 60*60, 'days': 60*60*24}

        return seconds[self.__time_unit]

    def total_seconds_per_day(self):

        return 86400
