import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Line plot for VM scaling
w_1 = np.array([104.511,
111.804,
114.441,
106.926,
110.788])


w_2 = np.array([119.307,
118.507,
118.148,
114.445])

w_3 = np.array([127.737,
121.178,
115.433,
121.989,
120.303])

w_4 = np.array([129.010,
129.862,
122.054,
125.859,
126.522])

w_2_4 = np.array([97.123,
97.628,
97.330,
93.627,
95.195])

w_4_4 = np.array([91.413,
91.008,
96.760,
90.725,
89.914])

w_2_2 = np.array([94.923,
91.293,
95.856,
97.769,
95.487])






# data = np.array([np.average(w_1),np.average(w_2),np.average(w_3),np.average(w_4)])
# print data
# plt.figure(1)
# x = range(1,5)
# plt.plot(x,data)
# plt.ylabel("Time taken for response (s)")

# plt.xlabel("# Workers")
# plt.xticks(x)
# plt.title("Time for response vs number of workers")
# plt.show()

# data = np.array([np.average(w_4),np.average(w_2_4),np.average(w_4_4)])
# print data
# plt.figure(3)
# left2 = []
# bar_width = 0.5
# space = 0.2
# colors = ['#6e95d3','#345c9b','#052454']
# labels = ['1 VM','2 VM', '4 VM']
# patch_list = []
# for i in range(len(colors)):
#     patch_list.append(mpatches.Patch(color=colors[i], label=labels[i]))
# bar_col = []
# for i in range(len(data)):
#     left2.append(bar_width*i + i/4*space)
#
# plt.bar(left2,data,width=bar_width, color=colors, edgecolor='k')
# plt.legend(handles=patch_list)
# plt.title("Difference in time whith 4 workers")
# plt.ylabel("Time for response (s)")
# plt.xlabel("Number of VMs")
# plt.xticks([])
# #plt.xticks([3.0/2*bar_width,3.0/2*bar_width+4*bar_width+space],["1 VM","2 VMs"])
# plt.show()

data = np.array([13.750,89.684,17.002,9.204])

data_lab = ['Problem 1 a I','Problem 1 b I', 'Problem 1 c I','Problem 1 b II']
bar_w = 0.6
color = '#0e8c10'
space = 0.3
left = []
lab_pos = []
for i in range(len(data)):
    left.append(i*(bar_w + space))
    lab_pos.append(i*(bar_w+space))
plt.figure(1)
plt.bar(left, data, width=bar_w, color = color)
plt.title("Computational time for individual problems")
plt.ylabel("Time taken (s)")
plt.xlabel("Problem")
plt.xticks(lab_pos,data_lab)

plt.show()






#plot for optimal
# plt.figure(4)
# data4 = np.array([60.725,28.858, 19.589, 15.842])
# x = range(1,5)
# ticklabs = ["1 VM", "2 VMs", "3 VMs", "4 VMs"]
# plt.plot(x,data4, color = '#5c037a')
# plt.title("Time for response when increasing number of VMs")
# plt.ylabel("Time for response")
# plt.xlabel("Number of VMs")
# plt.xticks(x,ticklabs)
# plt.show()



# 1 VM
# 1 worker
#
# [104.511,
# 111.804,
# 114.441,
# 106.926,
# 110.788]
#
# 2 worker
# [119.307,
# 118.507,
# 118.148,
# 114.445]
#
# 3 worker
# [127.737,
# 121.178,
# 115.433,
# 121.989,
# 120.303]
#
# 4 worker
# [129.010,
# 129.862,
# 122.054,
# 125.859,
# 126.522]
#
#
# 2 VM
# 1 worker
# [94.923,
# 91.293,
# 95.856,
# 97.769,
# 95.487]
#
# 2 worker
# [97.123,
# 97.628,
# 97.330,
# 93.627,
# 95.195]
#
# 4 VM
# 1 worker
# [91.413,
# 91.008,
# 96.760,
# 90.725,
# 89.914s]
