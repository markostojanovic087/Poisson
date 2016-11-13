# Poisson

Application solves the 3D Poisson Equation for n input sets of N x N x N size, where N has to be power of 2 and at least 32. For solving the equation we use FFT method. N is configured at build time, but n can be given in runtime as a command line argument. Input can be read from file or randomly generated. If user chooses so, all output and output convenient for plotting can be generated. The App runs both the DFE and the CPU implementation and compares respective times taken.
