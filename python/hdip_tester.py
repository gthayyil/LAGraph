import pygraphblas as gb
import laplacian as lap
from hdip_fiedler import hdip_fiedler
from operator import itemgetter
import time
"""
Hdip_tester: a program that is used to test hdip_fiedler, by loading in a matrix from ssget, modifying the matrix, applying hdip_fiedler to it, and using the result to partition the matrix.

Running this program will allow you to graph partition using hdip_fiedler on any matrix whose matrix id is stored in the matrixId list. 

Authors: Georgy Thayyil, Tim Davis
Acknowledgements: Michel Pelletier provided us with many helpful suggestions and assistance while developing this algorithm

The following cited websites were used to find a way to apply the permutation of the python sort.

Citations
    1:https://www.programiz.com/python-programming/methods/list/sort
    2:https://www.kite.com/python/answers/how-to-zip-two-lists-in-python
    3:https://www.kite.com/python/answers/how-to-unzip-a-list-of-tuples-in-python
"""
#matrixId is a list with all the matricies from the suitsparse collection that are symmetric, and only have 1 connected component.
matrixId=[    
    813, 1202, 1435, 2593, 1848, 409, 2639, 2597, 342, 362, 341, 1892, 785, 2623, 2585, 2629, 1296, 2604, 2769, 361, 2583, 2617, 953, 1590, 2596, 817, 1368, 2526, 786, 890, 2626, 2608,
    #2591, 1900, 340, 1452, 2518, 816, 1297, 788, 787, 789, 790, 1308, 2473, 2603, 855, 2461, 1231, 1383, 1919, 2605, 844, 845, 1251, 1434, 821, 822, 820, 1878, 352, 367, 1213, 2616, 804, 2778, 811, 2788, 1271,
    #1417, 1284, 2462, 1643, 54, 913, 52, 2619, 846, 2779, 2427, 964, 965, 966, 967, 2517, 536, 2595, 851, 1357, 1235, 2460, 1287, 928, 853, 2642, 1286, 2570, 852, 356, 761, 368, 1431, 1236, 1382, 2260, 1274,
    #973, 1282, 805, 2510, 849, 2474, 803, 854, 936, 2529, 1354, 1852, 2830, 1849, 1853, 1276, 2770, 1909, 1906, 1278, 1277, 2439, 850, 2283, 537, 856, 2373, 1208, 1209, 1450, 1275, 1283, 2624, 1898, 1910, 1423,
    2516, 1883, 1882, 1899, 1361, 1893, 1861, 1454, 1352, 2783, 1350, 2580, 2475, 2812, 859, 2519, 971, 937, 914, 858, 1421, 1349, 2813, 1355, 1356, 1290, 2265, 2266, 1351, 2571, 538, 1393, 1269, 2267, 2268, 
    1385, 1257, 1254, 1364, 2771, 1367, 1219, 1260, 369, 1264, 1285, 2831, 1265, 2476, 1857, 2572, 2483, 2514, 1258, 938, 860, 1362, 2498, 2513, 2486, 1580, 1581, 1582, 1583, 1584, 1585, 940, 941, 942, 943,
    944, 945, 946, 947, 948, 1353, 2487, 1398, 2579, 2488, 2329, 2577, 2578, 2512, 2581, 2477, 2385, 2509, 2573, 1411, 915, 2495, 1901, 939, 1856, 2464, 2463, 2387, 2772, 2458, 2484, 1267, 1455, 2485, 2386,
    2478, 1586, 2480, 2459, 2854, 2481, 2482, 2781, 1252, 2496, 2388, 1902, 916, 2773, 2479, 2511, 2544, 1903, 2782, 2774, 1904, 1905, 2775, 2776, 2780, 2856]

matrixId = [1883, 1882, 1899, 1361, 1893, 1861, 1454 ,1385 ,1257]
gb.options_set(burble=True)
#This for loop, loops through all the matrix ids in the list, and pauses after each one. 
for mid in matrixId:
    ht0=time.time()
    
    #The matrix is loaded in using ssget..
    print("--------------------------------------------------------------------------------------------------------------------")
    print("The matrix is loaded in using ssget..")
    print("--------------------------------------------------------------------------------------------------------------------")
    J=sorted(list(gb.Matrix.ssget(mid, binary_cache_dir="~/.ssgetpy")), key=itemgetter(0))[0]
    ht1=time.time()
    
    #The following lines of code set the matrix's type to FP32, makes sure the matrix is symmetric, and applies the laplacian to the modified matrix.
    print("--------------------------------------------------------------------------------------------------------------------")
    print("The following lines of code set the matrix's type to FP32, makes sure the matrix is symmetric, and applies the laplacian to the modified matrix.")
    print("--------------------------------------------------------------------------------------------------------------------")
    pattern=J[1].pattern(gb.types.FP32)
    h=pattern+pattern.transpose()
    omega=lap.laplacian(pattern)
    ht2=time.time()
    
    #Calls hdip_fiedler on the modified matrix
    print("--------------------------------------------------------------------------------------------------------------------")
    print("Calls hdip_fiedler on the modified matrix")
    print("--------------------------------------------------------------------------------------------------------------------")
    hdip2=hdip_fiedler(omega[0],omega[1])
    ht3=time.time()
    hdip_ht=ht3-ht2
    
    #assign hdip2[0] to a variable
    print("--------------------------------------------------------------------------------------------------------------------")
    print("assign hdip2[0] to a variable")
    print("--------------------------------------------------------------------------------------------------------------------")
    g=hdip2[0]
    h55=time.time()
    print(h55-ht3)

    #Converts g into a python list that transforms [1,2,3] into [(0,1),(1,2),(2,3)]
    print("--------------------------------------------------------------------------------------------------------------------")
    print("Converts g into a python list that transforms [1,2,3] into [(0,1),(1,2),(2,3)]")
    print("--------------------------------------------------------------------------------------------------------------------")
    z=list(g)
    h66=time.time()
    print(h66-h55)
    
    #Sort is applied with the itemgetter function to sort the ziped list based on the second value 
    print("--------------------------------------------------------------------------------------------------------------------")
    print("Sort is applied with the itemgetter function to sort the ziped list based on the second value")
    print("--------------------------------------------------------------------------------------------------------------------")
    z.sort(key=itemgetter(1))
    ht4=time.time()

    print("sort")
    print(ht4-h66)

    #The values of the indicies are extracted and stored into a list called j, then the lists length is taken, and the median index value is found.
    print("--------------------------------------------------------------------------------------------------------------------")
    print("The values of the indicies are extracted and stored into a list called j, then the lists length is taken, and the median index value is found.")
    print("--------------------------------------------------------------------------------------------------------------------")
    j=list(map(itemgetter(0), z))
    h5=time.time()
    print(h5-ht4)

    lenj=len(j)
    medind=int(lenj/2)
    h6=time.time()
    print(h6-h5)
    
    #A matrix based on the permutation of the applied sort to the list is extracted from the matrix that is loaded in by ssget.
    print("--------------------------------------------------------------------------------------------------------------------")
    print("A matrix based on the permutation of the applied sort to the list is extracted from the matrix that is loaded in by ssget.")
    print("--------------------------------------------------------------------------------------------------------------------")
    matrixS=J[1].extract_matrix(j,j)
    h7=time.time()
    print(h7-h6)
    
    #The matrix is partitioned, by using  matrix extract
    print("--------------------------------------------------------------------------------------------------------------------")
    print("The matrix is partitioned, by using  matrix extract")
    print("--------------------------------------------------------------------------------------------------------------------")
    matrixUpLef=matrixS.extract_matrix(slice(0,medind),slice(medind+1,lenj-1))
    h18=time.time()
    print(h18-h7)

    #The number of edgecuts are found by the number of elements in the partitioned matrix.
    print("--------------------------------------------------------------------------------------------------------------------")
    print("The number of edgecuts are found by the number of elements in the partitioned matrix.")
    print("--------------------------------------------------------------------------------------------------------------------")
    edgecut=matrixUpLef.nvals  

    #Used to record the total time it took to run everything except the following print statements
    ht8=time.time()
    print(ht8-h18)
    #The following are print statements which are used to compare the results of this code and hdip_fiedler with their counterparts.
    print("The hdip function result")
    print(hdip2)
    print("The first vector that is returned in the tuple")
    print("The edges cut is")
    print(edgecut)
    print("The value that is returned for lambda")
    print(hdip2[1])
    print("The following are the values inside the iters tuple")
    print(hdip2[2][0])
    print("The second value inside the iters tuple")
    print(hdip2[2][1]) 
    print("The id of the matrix")
    print(mid)
    print("The time of laplacian and other things before hdip except loading in time")
    print(ht2-ht1)
    print("The amount of time it takes to run Hdip_fiedler")
    print(hdip_ht)
    print("The time for zipper, meding etc, or everything after hdip")
    print(ht8-ht3)
    print("The time it takes to run my code excluding the time for the matrix to load in")
    print(ht8-ht1)
    print("The amount of time it takes to run everything except print statements, till the pause")
    print(ht8-ht0)
    print("(total time with ssget)/(time excluding ssget)")
    print((ht8-ht0)/(ht8-ht1))
    gb.options_get()
    wait = input("Press Enter to continue.")

