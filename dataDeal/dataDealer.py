from numpy import *
import random
#import numpy as np
#np.set_printoptions(suppress=True)
attributesCount = 22

#把标签与属性值分成D_setAll.txt，D_tagAll.txt两个文件
def classfiy(filename):
    readRow = []
    fs = open("D_setAll.txt","w")
    ft = open("D_tagAll.txt","w")
    for line in open(filename):
        readRow.clear()
        tmp = line.split('\t')
        for i in tmp:
            readRow.append(i.strip())
        for property in readRow[:-2]:
            fs.write(property+' ')
        fs.write(readRow[-2])
        fs.write('\n')
        ft.write(readRow[-1]+'\n')
        #print(readRow)
    fs.close()
    ft.close()

#标准化数据，生成数据在储存D_trainS.txt中
def standardization():
    readSet = []
    fs = open("D_setAll.txt")
    count = 0
    for line in fs.readlines():
        readRow = []
        count += 1
        tmp = line.split(' ')
        for i in tmp:
            if i == '':
                continue
            else:
                readRow.append(float(i.strip()))
        readSet.append(readRow)
    fs.close()
    dataSet = array(readSet)
    aver = array([0]*attributesCount)
    variance = array([0]*attributesCount)
    for i in range(count):
        aver = aver + dataSet[i]
    aver = aver/count
    print('平均值:',aver)
    tmpDataSet = zeros((count,attributesCount))
    finalResult = zeros((count,attributesCount))
    for i in range(count):
        tmpDataSet[i] = (dataSet[i] - aver)**2
    for i in range(count):
        variance = variance+tmpDataSet[i]
    variance = variance/count
    print('标准差:',variance)
    for i in range(count):
        finalResult[i] = (dataSet[i]-aver)/variance
    print(finalResult)
    fs2 = open("D_trainS.txt","w")
    for i in finalResult:
        for j in i[:-1]:
            fs2.write(str(j)+' ')
        fs2.write(str(i[-1]))
        fs2.write('\n')
    fs2.close()

#生成最后需要的八个文件，四个四个储存在original，standardization中
def getTrainTest():
    fr1 = open("D_trainS.txt")
    #fr2 = open("D_train_resultAll.txt")
    fr3 = open("D_setAll.txt")
    fr4 = open("D_tagAll.txt")
    datasetOAll =[]
    datasetAll = []
    datatagAll = []
    count = 0
    for line in fr1.readlines():
        readRow = []
        count += 1
        tmp = line.split(' ')
        for i in tmp:
            if i == '':
                continue
            else:
                readRow.append(float(i.strip()))
        datasetAll.append(readRow)

    for line in fr3.readlines():
        readRow =[]
        tmp = line.split(' ')
        for i in tmp:
            if i == '':
                continue
            else:
                readRow.append(float(i.strip()))
        datasetOAll.append(readRow)

    for line in fr4.readlines():
        readRow = []
        tmp = line.split(' ')
        for i in tmp:
            if i == '':
                continue
            else:
                readRow.append(float(i.strip()))
        datatagAll.append(readRow)
    print(datasetAll)
    print(datatagAll)
    randomList = random.sample(range(0,count),count)
    trainVolume = int(count*0.8)
    testVolume = count - trainVolume
    datasetOTrain = [0]*trainVolume
    datatagOTrain = [0]*trainVolume
    datasetOTest = [0]*testVolume
    datatagOTest = [0]*testVolume

    datasetTrain = [0]*trainVolume
    datatagTrain = [0]*trainVolume
    datasetTest = [0]*testVolume
    datatagTest = [0]*testVolume

    for i in range(trainVolume):
        datasetTrain[i] = datasetAll[randomList[i]]
        datatagTrain[i] = datatagAll[randomList[i]]
        datasetOTrain[i] = datasetOAll[randomList[i]]
        datatagOTrain[i] = datatagAll[randomList[i]]
    for i in range(testVolume):
        datasetTest[i] = datasetAll[randomList[i+trainVolume]]
        datatagTest[i] = datatagAll[randomList[i+trainVolume]]
        datasetOTest[i] = datasetOAll[randomList[i+trainVolume]]
        datatagOTest[i] = datatagAll[randomList[i+trainVolume]]
    fs10 = open("original/D_train.txt","w")
    fs20 = open("original/D_train_tag.txt","w")
    fs30 = open("original/D_test.txt","w")
    fs40 = open("original/D_test_tag.txt","w")
    for i in datasetOTrain:
        for j in i[:-1]:
            fs10.write(str(j)+' ')
        fs10.write(str(i[-1]))
        fs10.write('\n')
    fs10.close()
    for i in datatagOTrain:
        #for j in i[:-1]:
        #    fs2.write(str(j)+' ')
        fs20.write(str(i[0]))
        fs20.write('\n')
    fs20.close()
    for i in datasetOTest:
        for j in i[:-1]:
            fs30.write(str(j)+' ')
        fs30.write(str(i[-1]))
        fs30.write('\n')
    fs30.close()
    for i in datatagOTest:
        #for j in i[:-1]:
        #    fs4.write(str(j)+' ')
        fs40.write(str(i[0]))
        fs40.write('\n')
    fs40.close()
    fs1 = open("standardization/D_train.txt","w")
    fs2 = open("standardization/D_train_tag.txt","w")
    fs3 = open("standardization/D_test.txt","w")
    fs4 = open("standardization/D_test_tag.txt","w")

    for i in datasetTrain:
        for j in i[:-1]:
            fs1.write(str(j)+' ')
        fs1.write(str(i[-1]))
        fs1.write('\n')
    fs1.close()
    for i in datatagTrain:
        #for j in i[:-1]:
        #    fs2.write(str(j)+' ')
        fs2.write(str(i[0]))
        fs2.write('\n')
    fs2.close()
    for i in datasetTest:
        for j in i[:-1]:
            fs3.write(str(j)+' ')
        fs3.write(str(i[-1]))
        fs3.write('\n')
    fs3.close()
    for i in datatagTest:
        #for j in i[:-1]:
        #    fs4.write(str(j)+' ')
        fs4.write(str(i[0]))
        fs4.write('\n')
    fs4.close()
if __name__ == "__main__":
    classfiy("mushroom.txt")
    standardization()
    getTrainTest()


