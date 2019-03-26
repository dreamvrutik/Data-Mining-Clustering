import random

#function to parse the dataset
def read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

#Distance function
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

headers = [] #stores the headers of amino acids
sequences = [] #stores the sequences

with open('aminoacids.fa') as fp:
    for name, seq in read_fasta(fp):
        headers.append(name)
        sequences.append(seq)

distance = []

for i in range(0,len(sequences)):
    temp = []
    for j in range(0,len(sequences)):
        temp.append([])
    distance.append(temp)

#Compute distance between every pair of points
for i in range(0,len(sequences)):
    for j in range(0,len(sequences)):
        print(i,j)
        distance[i][j]=editDistance(sequences[i],sequences[j],len(sequences[i]),len(sequences[j]))


points = [x for x in range(0,311)]

k = 6
centroids = [] #stores the k centroids
clusters = [] #stores the points in the k clusters

#Randomly selecting initial k centroids
for i in range(0,k):
    x = random.choice(points)
    points.remove(x)
    centroids.append(x)


print(centroids)
#K-means algorithm
for run in range(0,10):
    clusters = []
    for i in range(0,k):
        clusters.append([])
    #forming clusters  
    for i in range(0,len(sequences)):
        min_dist = 1000000000000000
        min_cluster = -1
        for j in range(0,k):
            temp_dist = distance[i][centroids[j]]
            if(temp_dist < min_dist):
                min_dist = temp_dist
                min_cluster = j
        clusters[min_cluster].append(i)
        
    flag = 1
    
    #computing new centroids
    for i in range(0,k):
        min_sum = 1000000000000000
        new_centroid = -1
        for x in clusters[i]:
            dist_sum = 0
            for y in clusters[i]:
                dist_sum+=distance[x][y]
            #print(x,dist_sum)
            if dist_sum < min_sum:
                min_sum = dist_sum
                new_centroid = x
        if new_centroid != centroids[i]:
            flag = 0
        centroids[i] = new_centroid
    
    print(centroids)
    if(flag):
        break
print(clusters)
print(centroids)
    

