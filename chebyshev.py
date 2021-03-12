#!/usr/bin/env python
# coding: utf-8

# In[32]:


import math

pi = math.pi

def studentInfo():
    '''prints student information
    '''
    print("""
        NAME: ADERONMU AYOMIDE M.
        MATRIC NO: 160403029
        DEPARTMENT: ELECTRICAL/ELECTRONIC ENGINEERING
        COURSE CODE: EEG 416
        COURSE TITLE: ACTIVE NETWORKS AND SYNTHESIS
        PROGRAM: GENERATION OF COEFFICIENTS OF CHEBYSHEV APPROXIMATION FUNCTIONS
        CHEBYSHEV APPROXIMATION TABLE FOR FILTER NETWORKS
        PACKAGE USED: Python
    """, '\n')


def calc_epsilon(Apmax):
    e = math.sqrt(10**(0.1*Apmax)-1)
    return e


def calc_order(Apmax, Asmin, wp, ws):
    """
    Parameters:
    Apmax - maximum passand loss
    Asmin - minimum stopband loss
    wp - passband frequency
    ws - stopband frequency
    
    ----
    output: order of the required chebyshev filter
    """

    n = 0
    A = calc_epsilon(Asmin) / calc_epsilon(Apmax)
    omega = ws/wp

    if omega > 1:
        n = math.acosh(A)/math.acosh(omega)
    else:
        n = math.acos(A)/math.acos(omega)

    n = math.ceil(n)
    return n



def calculate_rho_k(n, k, eps):
    """
    Parameters:
    n: order
    k: iteration
    eps: epselom
    
    ---
    output: rho
    """
    rho = math.sin((pi/2) * (1 + 2 * k) / n) * math.sinh(math.asinh(1/eps) / n)
    
    return rho


def calculate_w_k(n, k, eps):
    """
    Parameters:
    n: order
    k: iteration
    eps: epselom
    
    ---
    output: omega
    """
    omega = math.cos((pi/2) * (1 + 2 * k) / n)  *  math.cosh(math.asinh(1/eps) / n)
    
    return omega


def calc_quadratic(rho, w):
    """
    Parameters:
    roots of |H(s)|-square: rho and omega
        
    ---
    output - quadratic polynomial part of the Numerator of H(s)
    """
    M = round(2*rho, 5)
    N = round(rho*rho + w*w, 5)

    return "(s^2 + {}s + {})".format(M, N)


def find_coeffs(Apmax, n):
    """
    Apmax - maximum passand loss
    Asmin: minimum stopband loss
    wp - passband frequency
    ws - stopband frequency
    
    ---
    output: a string of the numerator, denominator and constant of the chebyshev filter
    """

    epsilon = calc_epsilon(Apmax)
    denominator = 1             # initalize the denominator value
    numerator = ''

    # calculate the resulting quadratic while within the range of order, n//2
    for k in range((n//2)):
        rho = calculate_rho_k(n, k, epsilon)
        w = calculate_w_k(n, k, epsilon)
        denominator *= (rho*rho + w*w)
        numerator += calc_quadratic(rho,w)

    # calculate first degree polynomial of the Numerator
    if(n%2 == 1):
        k = n//2
        rho = round(calculate_rho_k(n, k, epsilon), 5)
        denominator *= rho
        w = calculate_w_k(n, k, epsilon)
        numerator += "(s + {})".format(rho)

    else:
        denominator *= 10**(-1*Apmax/20)

    denominator = round(denominator, 5)
    print("order = ", n, "\nNumerator of H(s)= ", numerator, "\nDenominator Constant, K = ",denominator,"\n")
    return "{},{},{}\n".format(n, numerator, denominator)




# PRINTING OUT THE RESULT and WRITING IT INTO AN CSV FILE

with open("chebyshev.csv","w") as file:
    
    studentInfo()       # print student info
    Apmax = [0.25, 0.50, 1.00]
    for ap in Apmax:
        print(f"__________________________________CHEBYSHEV VALUES FOR Apmax={ap} dB_________________________________________")
        file.write(f"CHEBYSHEV VALUES FOR Apmax={ap} dB\n")
        file.write("n,Numerator of H(s),Denominator Constant K\n")
        for i in range(1,11):
            file.write(find_coeffs(ap,i))

    file.write("\n\n\n\n")
    print("\n\n\n\n")
    
    

# Osemudiamen Itua (https://github.com/Ose-4g) contributed heavily to the development of this program

