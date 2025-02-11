#!/usr/bin/env python

BRIGHT_BLACK = '\033[90m'
BRIGHT_RED = '\033[91m'
BRIGHT_GREEN = '\033[92m'
BRIGHT_YELLOW = '\033[93m'
BRIGHT_BLUE = '\033[94m'
BRIGHT_MAGENTA = '\033[95m'
BRIGHT_CYAN = '\033[96m'
BRIGHT_WHITE = '\033[97m'
BRIGHT_END = '\033[0m'


class C:
    K = '\033[30m'  #BLACK
    R = '\033[31m'  #RED
    G = '\033[32m'  #GREEN
    Y = '\033[33m'  #YELLOW
    B = '\033[34m'  #BLUE 
    M = '\033[35m'  #MAGENTA
    C = '\033[36m'  #CYAN 
    W = '\033[37m'  #WHITE
    U = '\033[4m'   #UNDERLINE
    Z = '\033[0m'   #END(RESET)

if __name__ == '__main__':
    print(C.R + 'TEST' + C.Z)
    print(BRIGHT_YELLOW + 'TEST' + BRIGHT_END)
    print(C.R + 'TE' + C.B + 'ST' + C.Z)
    print(C.U + 'TEST' + C.Z + 'ST' + C.Z)
