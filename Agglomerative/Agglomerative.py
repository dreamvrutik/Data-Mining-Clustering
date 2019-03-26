# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 03:09:34 2019

@author: Kunal
"""

#Library to plot graphs

import matplotlib.pyplot as plt


#Edit distance between strings

def editDistance(str1, str2, m, n): 
    # Create a table to store results of subproblems 
    dp = [[0 for x in range(n+1)] for x in range(m+1)] 
  
    # Fill d[][] in bottom up manner 
    for i in range(m+1): 
        for j in range(n+1): 
  
            # If first string is empty, only option is to 
            # insert all characters of second string 
            if i == 0: 
                dp[i][j] = j    # Min. operations = j 
  
            # If second string is empty, only option is to 
            # remove all characters of second string 
            elif j == 0: 
                dp[i][j] = i    # Min. operations = i 
  
            # If last characters are same, ignore last char 
            # and recur for remaining string 
            elif str1[i-1] == str2[j-1]: 
                dp[i][j] = dp[i-1][j-1] 
  
            # If last character are different, consider all 
            # possibilities and find minimum 
            else: 
                dp[i][j] = 1 + min(dp[i][j-1],        # Insert 
                                   dp[i-1][j],        # Remove 
                                   dp[i-1][j-1])    # Replace 
  
    return dp[m][n] 




#Importing the amino acid sequences

sequences = []
with open('AminoAcids.txt') as fp:
    contents = fp.readlines()
    current_sequence = ""
    count = 0
    for line in contents:
        if line[0] == '>':
           count = count + 1
           if count>1:
               sequences.append(current_sequence[:-1])
           current_sequence = ""
        else:
            current_sequence = current_sequence + line[:-1]
fp.close()
    



#Creating Initial Distance Matrix
          
n = len(sequences)
distance_matrix = [[0 for x in range(0,n)] for y in range(0,n)]
for i in range(0,n):
    for j in range(0,n):
        distance_matrix[i][j] = editDistance(sequences[i],sequences[j],len(sequences[i]),len(sequences[j]))

        




#Single Linkage Clusters 
     
temp = [[0 for x in range(0,n)] for y in range(0,n)]
costs = []
x_axis = []
y_axis = []
plt.title("Single Linkage Agglomerative Clustering Of Amino Acid Sequences")
plt.xlabel("Cluster Index")
plt.ylabel("Distance")
for i in range(0,n):
    x_axis.append(i+1)
    y_axis.append(0)
for i in range(0,n):
    for j in range(0,n):
        temp[i][j] = distance_matrix[i][j]
linkages = []
for i in range(0,n-1):
    val = 1000000000
    posx = n
    posy = n
    for j in range(0,n):
        for k in range(0,n):
            if k == j:
                continue
            if temp[j][k] == -1:
                continue
            if temp[j][k] < val:
                val = temp[j][k]
                posx = j
                posy = k
    if posx > posy:
        posx,posy = posy,posx
    X = [x_axis[posx],x_axis[posx]]
    Y = [y_axis[posx],val]
    plt.plot(X,Y)
    X = [x_axis[posy],x_axis[posy]]
    Y = [y_axis[posy],val]
    plt.plot(X,Y)
    y_axis[posx]=val
    y_axis[posy]=val
    X = [x_axis[posx],x_axis[posy]]
    Y = [val,val]
    plt.plot(X,Y)
    x_axis[posx] = (x_axis[posx]+x_axis[posy])/2
    tup = (posx,posy)
    linkages.append(tup)
    costs.append(val)
    for j in range(0,n-1):
        val = min(temp[posx][j],temp[posy][j])
        temp[posx][j]=val
        temp[j][posx]=val
    for j in range(0,n-1):
        temp[posy][j] = -1
        temp[j][posy] = -1
plt.show()





#Complete Linkage Clusters
     
costs = []
x_axis = []
y_axis = []
plt.title("Complete Linkage Agglomerative Clustering Of Amino Acid Sequences")
plt.xlabel("Cluster Index")
plt.ylabel("Distance")
for i in range(0,n):
    x_axis.append(i+1)
    y_axis.append(0)
for i in range(0,n):
    for j in range(0,n):
        temp[i][j] = distance_matrix[i][j]
linkages_complete = []
for i in range(0,n-1):
    val = 1000000000
    posx = n
    posy = n
    for j in range(0,n):
        for k in range(0,n):
            if k == j:
                continue
            if temp[j][k] == -1:
                continue
            if temp[j][k] < val:
                val = temp[j][k]
                posx = j
                posy = k
    if posx > posy:
        posx,posy = posy,posx
    X = [x_axis[posx],x_axis[posx]]
    Y = [y_axis[posx],val]
    plt.plot(X,Y)
    X = [x_axis[posy],x_axis[posy]]
    Y = [y_axis[posy],val]
    plt.plot(X,Y)
    y_axis[posx]=val
    y_axis[posy]=val
    X = [x_axis[posx],x_axis[posy]]
    Y = [val,val]
    plt.plot(X,Y)
    x_axis[posx] = (x_axis[posx]+x_axis[posy])/2
    tup = (posx,posy)
    linkages_complete.append(tup)
    costs.append(val)
    for j in range(0,n-1):
        val = max(temp[posx][j],temp[posy][j])
        temp[posx][j]=val
        temp[j][posx]=val
    for j in range(0,n-1):
        temp[posy][j] = -1
        temp[j][posy] = -1
plt.show()





#Group Average Clusters
     
costs = []
x_axis = []
y_axis = []
cluster_size = []
plt.title("Group Average Agglomerative Clustering Of Amino Acid Sequences")
plt.xlabel("Cluster Index")
plt.ylabel("Distance")
for i in range(0,n):
    x_axis.append(i+1)
    y_axis.append(0)
    cluster_size.append(1)
for i in range(0,n):
    for j in range(0,n):
        temp[i][j] = distance_matrix[i][j]
linkages_average = []
for i in range(0,n-1):
    val = 1000000000
    posx = n
    posy = n
    for j in range(0,n):
        for k in range(0,n):
            if k == j:
                continue
            if temp[j][k] == -1:
                continue
            if temp[j][k] < val:
                val = temp[j][k]
                posx = j
                posy = k
    if posx > posy:
        posx,posy = posy,posx
    X = [x_axis[posx],x_axis[posx]]
    Y = [y_axis[posx],val]
    plt.plot(X,Y)
    X = [x_axis[posy],x_axis[posy]]
    Y = [y_axis[posy],val]
    plt.plot(X,Y)
    y_axis[posx]=val
    y_axis[posy]=val
    X = [x_axis[posx],x_axis[posy]]
    Y = [val,val]
    plt.plot(X,Y)
    x_axis[posx] = (x_axis[posx]+x_axis[posy])/2
    tup = (posx,posy)
    linkages_average.append(tup)
    costs.append(val)
    for j in range(0,n-1):
        if(temp[posx][j] == -1):
            continue
        val = (temp[posx][j]*cluster_size[posx]+temp[posy][j]*cluster_size[posy])/((cluster_size[posx]+cluster_size[posy])*cluster_size[j])
        temp[posx][j]=val
        temp[j][posx]=val
    cluster_size[posx] = cluster_size[posx] + cluster_size[posy]
    cluster_size[posy] = 0
    for j in range(0,n-1):
        temp[posy][j] = -1
        temp[j][posy] = -1
plt.show()