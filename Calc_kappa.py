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

def targetFunctionSigmaZ(theta, phi):
    sint = math.sin(theta)
    return pow(sint,3)

def targetFunctionSigmaX(theta, phi):
    sint = math.sin(theta)
    cosp = math.cos(phi)
    sint2 = pow(sint,2)
    cosp2 = pow(cosp,2)
    return pow(1.0 - sint2*cosp2, 1.5)

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

def integrationSigmaZ(k):
    nTheta = 30
    nPhi = 30
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
            f1 = targetFunctionSigmaZ(theta1,phi1)
            f2 = targetFunctionSigmaZ(theta1,phi2)
            sum1 += 0.5 * dPhi * (f1 + f2) * pdfXSinTheta(theta1, k)
            f1 = targetFunctionSigmaZ(theta2,phi1)
            f2 = targetFunctionSigmaZ(theta2,phi2)
            sum2 += 0.5 * dPhi * (f1 + f2) * pdfXSinTheta(theta2, k)
        result += 0.5 * dTheta * (sum1 + sum2)
    return result

def integrationSigmaX(k):
    nTheta = 30
    nPhi = 30
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
            f1 = targetFunctionSigmaX(theta1,phi1)
            f2 = targetFunctionSigmaX(theta1,phi2)
            sum1 += 0.5 * dPhi * (f1 + f2) * pdfXSinTheta(theta1, k)
            f1 = targetFunctionSigmaX(theta2,phi1)
            f2 = targetFunctionSigmaX(theta2,phi2)
            sum2 += 0.5 * dPhi * (f1 + f2) * pdfXSinTheta(theta2, k)
        result += 0.5 * dTheta * (sum1 + sum2)
    return result

def f(k, ans):
    return integrationSigmaX(k)/integrationSigmaZ(k) - ans

sigma_x = 0.0038610038610 # conductivty (S/m) along x-axis
sigma_z = 0.0002201673272 # conductivty (S/m) along z-axis
sigma_l = 21.0 # conductivity (S/m) of the coductive phase in cracks

ans = sigma_x / sigma_z
eps = 0.001
a = -10.0
b = 10.0
n = 1
while True:
    c = (a + b)/2
    if f(a, ans) * f(c, ans) < 0:
        b = c
    else:
        a = c
    if abs(f(c, ans)) < eps:
        fileName = 'kappa_phi.txt'
        f = open(fileName, 'w')
        f.write("kappa: " + str(c) + "\n")
        f.write("crack volume fraction: " + str(sigma_z/integrationSigmaZ(c)/sigma_l) + "\n")
        f.close()
        break
    n += 1
