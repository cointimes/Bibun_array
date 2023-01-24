#changed from cube_analysis_class.py
#this program may cost 30sec.
import numpy as np
import matplotlib.pyplot as plt
#import time
#import re, os, sys


def sabun_a(n_kai, array2dim):  #return(array_xyz) 3dim
    """calculation of the derivation, please use under some functions"""
    cell = np.array([[2*np.pi/200, 0.], [0., 2*np.pi/200]]) #
    [dx, dy] = np.power(np.sqrt(np.sum(np.square(cell), axis=1)), n_kai)
    a_size = array2dim.shape
    def weight_p(num_o, num_p):
        """num_order and num_point is int"""
        d3  = np.array([-1., 0., 1.])/2.
        d5  = np.array([1., -8., 0., 8., -1.])/12.
        d7  = np.array([-1., 9., -45., 0., 45., -9., 1.])/60.
        d9  = np.array([3., -32., 168., -672., 0., 672., -168., 32., -3.])/840.
        dd3 = np.array([1., -2., 1.])
        dd5 = np.array([-1., 16., -30., 16., -1.])/12.
        dd7 = np.array([2., -27., 270., -490., 270., -27., 2.])/180.
        dd9 = np.array([-9., 128., -1008., 8064., -14350., 8064., -1008., 128., -9.])/5040.
        weight_list = [[],[d3,d5,d7,d9],[dd3,dd5,dd7,dd9]]
        return(weight_list[num_o][num_p//2 -1])
    NUM_POINT = 9
    a_WEIGHT = weight_p(n_kai,NUM_POINT)        # select coefficients
    a_out = []                              # finally change shape to (i,j,k,3)
    for i in range(array2dim.shape[0]):
        for j in range(array2dim.shape[1]):
            #for k in range(array_xyz.shape[2]):
            #from ijk=000 to ijk=199,199,199
                #ls_p = [0,0,-2:3] ~ [199,199,197:2]
                #ls_kp = array_xyz[i,j,k-NUM_POINT//2:k+1+NUM_POINT//2]
                ls_jp = array2dim[i,j-NUM_POINT//2:j+1+NUM_POINT//2]
                ls_ip = array2dim[i-NUM_POINT//2:i+1+NUM_POINT//2,j]
                # exception handling
                # map is not much faster
                #if (len(ls_kp) < NUM_POINT) :
                #    ls_kp = np.array(list(map(lambda p: array2dim[i,j,p], list(map(lambda q: q%a_size[2], np.arange(k-NUM_POINT//2, k+1+NUM_POINT//2))))))
                if (len(ls_jp) < NUM_POINT) :
                    ls_jp = np.array(list(map(lambda p: array2dim[i,p], list(map(lambda q: q%a_size[1], np.arange(j-NUM_POINT//2, j+1+NUM_POINT//2))))))
                if (len(ls_ip) < NUM_POINT) :
                    ls_ip = np.array(list(map(lambda p: array2dim[p,j], list(map(lambda q: q%a_size[0], np.arange(i-NUM_POINT//2, i+1+NUM_POINT//2))))))
                a_out.extend([np.dot(ls_ip, a_WEIGHT)/dx, np.dot(ls_jp, a_WEIGHT)/dy])
    a_out = np.array(a_out).reshape(array2dim.shape + (2,))
    a_out1 = np.sum(a_out, axis=2)
    return(a_out1)



def imaging_sabun(array):
    plotarray = [array, sabun_a(1,array), sabun_a(2,array)]
    tit = ['f(x,y)', 'd/dx + d/dy: meeningless', 'd/dx^2 + d/dy^2']
    fig = plt.figure(figsize=(12., 3.0), constrained_layout=True)
    axs=fig.subplots(1,3, squeeze=False)
    axs = axs.flatten()
    for i in range(3):
        mm=axs[i].imshow(plotarray[i], origin='lower')
        axs[i].set_title(tit[i])
        fig.colorbar(mm)
    plt.pause(5)
    #plt.savefig('bibun_imaging.png', format='png')
    plt.close()
    return('ok')



array22_0 = np.zeros((200,200), dtype=np.float32)
# 0 <= x < 200, 0 <= y < 200
# z[x, y] = 1
for i in range(200):
    for j in range(200):
        array22_0[i,j] = 1.

array22_1 = np.zeros((200,200), dtype=np.float32)
# z[x, x^2 <= y < x^2+c] = 1
for i in range(200):
    for j in range(200):
        if(     i <=(j-100)**2 /50 +10
            and i > (j-100)**2 /50 ):
            array22_1[i,j] = 1.

array22_2 = np.zeros((200,200), dtype=np.float32)
# z[x, sin(x) <= y < sin(x)+c ]
for i in range(200):
    for j in range(200):
        if(     i <= np.sin(j*2*np.pi/200)*100+105
            and i >  np.sin(j*2*np.pi/200)*100+ 95):
            array22_2[i,j] = 1.

array22_3 = np.zeros((200,200), dtype=np.float32)
# z[x, ]
for i in range(200):
    for j in range(200):
        #print(i, j, array22_3[i,j])
        array22_3[i,j] = np.exp(np.abs((i-100)*2*np.pi/200)+np.abs((j-100)*2*np.pi/200))
        #print(array22_3[i,j])  -> inf

array22_4 = np.zeros((200,200), dtype=np.float32)
for i in range(200):
    for j in range(200):
        array22_4[i,j] = np.sin(i*2*np.pi/200)

array22_5 = np.zeros((200,200), dtype=np.float32)
for i in range(200):
    for j in range(200):
        array22_5[i,j] = np.sin((i-100)*2*np.pi/200)*np.sin((j-100)*2*np.pi/200)


def main():
    imaging_sabun(array22_0)
    #imaging_sabun(array22_1)
    imaging_sabun(array22_2)
    #imaging_sabun(array22_3)
    imaging_sabun(array22_4)
    imaging_sabun(array22_5)

if __name__ == '__main__':
    main()