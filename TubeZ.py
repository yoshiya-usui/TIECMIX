# MIT License
#
# Copyright (c) 2024 Yoshiya Usui
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import math
import scipy.special as spc

def targetFunction(theta, phi):
    cost = math.cos(theta)
    return pow(cost, 3)

def pdf(theta, k):
    cost = math.cos(theta)
    sint = math.sin(theta)
    return math.exp( k * pow(cost,2) ) * 2.0 / ( math.pi * spc.hyp1f1(0.5, 1.5, k) )

def pdfXSinTheta(theta, k):
    return pdf(theta, k) * math.sin(theta)

def check(nTheta, nPhi, k):
    dTheta = 0.5 * math.pi / float(nTheta)
    dPhi = 0.5 * math.pi / float(nPhi)
    result = 0.0
    for iTheta in range(nTheta):
        theta1 = dTheta * float(iTheta)
        theta2 = dTheta * float(iTheta+1)
        f1 = pdfXSinTheta(theta1, k)
        f2 = pdfXSinTheta(theta2, k)
        result += 0.5 * dTheta * (f1 + f2)
    result *= 0.5 * math.pi
    return result

def integration(k):
    nTheta = 200
    nPhi = 200
    # print(check(nTheta, nPhi, k))
    dTheta = 0.5 * math.pi / float(nTheta)
    dPhi = 0.5 * math.pi / float(nPhi)
    result = 0.0
    for iTheta in range(nTheta):
        theta1 = dTheta * float(iTheta)
        theta2 = dTheta * float(iTheta+1)
        sum1 = 0.0
        sum2 = 0.0
        for iPhi in range(nPhi):
            phi1 = dPhi * float(iPhi)
            phi2 = dPhi * float(iPhi+1)
            f1 = targetFunction(theta1,phi1)
            f2 = targetFunction(theta1,phi2)
            sum1 += 0.5 * dPhi * (f1 + f2) * pdfXSinTheta(theta1, k)
            f1 = targetFunction(theta2,phi1)
            f2 = targetFunction(theta2,phi2)
            sum2 += 0.5 * dPhi * (f1 + f2) * pdfXSinTheta(theta2, k)
        result += 0.5 * dTheta * (sum1 + sum2)
    return result

phi=0.01 #volume fraction (0 <= phi <= 1.0)
sigma_l = 1.0 # conductivity (S/m) of the coductive phase in cracks
fileName = 'tube_sigma_z.csv'
f = open(fileName, 'w')
f.write("kappa, sigma_z"+'\n')
list = (100.0, 50.0, 30.0, 20.0, 10.0, 5.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -5.0, -10.0, -20.0, -30.0, -50.0, -100.0)
for kappa in list:
    sigma_z = integration(kappa) * sigma_l * phi # conductivity (S/m)
    s = str(kappa)+','+str(sigma_z)+'\n'
    f.write(s)
    f.flush()
f.close()
