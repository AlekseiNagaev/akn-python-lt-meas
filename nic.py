import PyDAQmx as nidaq
import numpy as np
import time
import matplotlib.pyplot as plt
t = nidaq.Task()
t.CreateAIVoltageChan("Dev1/ai0", None, nidaq.DAQmx_Val_Diff, 0, 0.1, nidaq.DAQmx_Val_Volts, None)
t.CfgSampClkTiming("", int(1e6), nidaq.DAQmx_Val_Rising, nidaq.DAQmx_Val_FiniteSamps, int(20e6))
t.StartTask()
time.sleep(20)
data = np.zeros((int(20e6),), dtype=np.float64)
read = nidaq.int32()
t.ReadAnalogF64(int(20e6), 25, nidaq.DAQmx_Val_GroupByChannel, data, len(data), nidaq.byref(read), None)
t.StopTask()
