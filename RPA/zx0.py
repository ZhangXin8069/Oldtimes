# count = 0
#
#
# # 㐀-鿕
#
#
# def function3():
#     import numpy as np
#     import sklearn.cluster as skc
#     from sklearn import metrics
#     import matplotlib.pyplot as plt
#
#     def f(x):
#         mac2id = dict()
#         onlinetimes = []
#         f = open(
#             r'D:\Programming\Python\txt&xlsx&docx\TestData.txt',
#             encoding='utf-8')
#         for line in f:
#             mac = line.split(',')[2]
#             onlinetime = int(line.split(',')[6])
#             starttime = int(line.split(',')[4].split(' ')[1].split(':')[0])
#             if mac not in mac2id:
#                 mac2id[mac] = len(onlinetimes)
#                 onlinetimes.append((starttime, onlinetime))
#             else:
#                 onlinetimes[mac2id[mac]] = [(starttime, onlinetime)]
#         real_X = np.array(onlinetimes).reshape((-1, 2))
#
#         X = real_X[:, 0:1]
#
#         db = skc.DBSCAN(eps=0.01, min_samples=20).fit(X)
#         labels = db.labels_
#
#         print('Labels:')
#         print(labels)
#         raito = len(labels[labels[:] == -1]) / len(labels)
#         print('Noise raito:', format(raito, '.2%'))
#
#         n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
#
#         print('Estimated number of clusters: %d' % n_clusters_)
#         print("Silhouette Coefficient: %0.3f" %
#               metrics.silhouette_score(X, labels))
#
#         for i in range(n_clusters_):
#             print('Cluster ', i, ':')
#             print(list(X[labels == i].flatten()))
#
#         plt.hist(X, 24)
#
#         pass
#
#
# def function2():
#     import numpy as np
#     from sklearn.cluster import KMeans
#
#     def loadData(filePath):
#         fr = open(filePath, 'r+', encoding='utf-8')
#         lines = fr.readlines()
#         retData = []
#         retCityName = []
#         for line in lines:
#             items = line.strip().split(",")
#             retCityName.append(items[0])
#             retData.append([float(items[i]) for i in range(1, len(items))])
#         return retData, retCityName
#
#     if __name__ == '__main__':
#         data, cityName = loadData(
#             r'D:\Programming\Python\txt&xlsx&docx\city.txt')
#         km = KMeans(n_clusters=4)
#         label = km.fit_predict(data)
#         expenses = np.sum(km.cluster_centers_, axis=1)
#         # print(expenses)
#         CityCluster = [[], [], [], []]
#         for i in range(len(cityName)):
#             CityCluster[label[i]].append(cityName[i])
#         for i in range(len(CityCluster)):
#             print("Expenses:%.2f" % expenses[i])
#             print(CityCluster[i])
#     pass
#
#
# def function1():
#     import numpy as np
#     import sklearn.cluster as skc
#     from sklearn import metrics
#     import matplotlib.pyplot as plt
#     mac2id = dict()
#     onlinetimes = []
#     f = open(
#         r'D:\Programming\Python\txt&xlsx&docx\TestData.txt',
#         encoding='utf-8')
#
#     for line in f:
#
#         # '''读取每条数据中的mac地址，开始上网时间，上网时长'''
#         mac = line.split(',')[2]
#         onlinetime = int(line.split(',')[6])
#         starttime = int(line.split(',')[4].split(' ')[1].split(':')[0])
#
#         # '''mac2id是一个字典，key是mac地址，value是对应mac地址的上网时长以及开始上网时间'''
#         if mac not in mac2id:
#             mac2id[mac] = len(onlinetimes)
#             onlinetimes.append((starttime, onlinetime))
#         else:
#             onlinetimes[mac2id[mac]] = [(starttime, onlinetime)]
#     real_X = np.array(onlinetimes).reshape((-1, 2))
#     # '''调用DBSCAN方法进行训练，labels存放每个数据对应的簇标签'''
#     X = real_X[:, 0:1]
#     db = skc.DBSCAN(eps=0.01, min_samples=20).fit(X)
#     labels = db.labels_
#     # '''打印数据被记上的标签，计算标签为-1的数据个数，得出噪声数据的比例'''
#     print('Labels:')
#     print(labels)
#     raito = len(labels[labels[:] == -1]) / len(labels)
#     print('Noise raito:', format(raito, '.2%'))
#     # '''计算簇的个数并打印，评价聚类效果'''
#     n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
#     print('Estimated number of clusters: %d' % n_clusters_)
#     print("Silhouette Coefficient: %0.3f" %
#           metrics.silhouette_score(X, labels))
#
#     # '''打印各簇标号以及各簇内数据'''
#
#     for i in range(n_clusters_):
#         print('Cluster ', i, ':')
#         print(list(X[labels == i].flatten()))
#     plt.hist(X, 24)
#
#     plt.xlabel("time", fontsize=24)
#     plt.ylabel("student's number", fontsize=24)
#     plt.show()
#     pass
#
#
# def function0():
#     import numpy as np
#     from sklearn.cluster import KMeans
#
#     def loadData(filePath):
#         # 打开文件流
#         fr = open(filePath, 'r+', encoding='UTF-8')
#         # 完整的读取，一次性读取整个文件，按行读取
#         lines = fr.readlines()
#         # 将读入的数据进行拆分，分为数据和城市名
#         retData = []
#         retCityName = []
#         for line in lines:
#             # 去除字符串首尾的空格或者回车，并使用“，”进行分割
#             items = line.strip().split(",")
#             # 每一行的开头是城市名臣
#             retCityName.append(items[0])
#             # 将数据组合成一个列表，并且强制转换类型为float浮点型
#             retData.append([float(items[i]) for i in range(1, len(items))])
#         return retData, retCityName
#
#     if __name__ == '__main__':
#         # 使用读取数据，获取城市名和相关的数据
#         data, cityName = loadData(
#             'D:/Programming/Python/txt&xlsx&docx/city.txt')
#         # 创建指定簇数量KMeans对象实例
#         km = KMeans(n_clusters=6)
#         # 加载数据，进行训练，获得标签，总共是四个簇，就是四个标签，将给31个数据，每个数据都打上0-3的标签
#         label = km.fit_predict(data)
#         # 计算出每一个簇形成的所有的行内的数据，计算出该簇内的数据的和
#         expenses = np.sum(km.cluster_centers_, axis=1)
#         # 总共是四个标签，四个集合，按照打上的标签将城市名进行分类
#         CityCluster = [[], [], [], [], [], []]
#         # 遍历所有的标签，并将对应的城市根据标签加上对应的簇中
#         for i in range(len(cityName)):
#             CityCluster[label[i]].append(cityName[i])
#         # 遍历所有的簇中心的数量，总共就只有四个，进行打印输出
#         for i in range(len(CityCluster)):
#             print("Expenses:%.2f" % expenses[i])
#             print(CityCluster[i])
#
#
# def _每日新闻():
#     import requests
#     import re
#     url = 'https://news.china.com/international/'
#     r = requests.get(url)
#     r.encoding = 'utf-8'
#     html = r.text
#     alist = re.findall(
#         r'.html" target="_blank">.+</a></h3><span class="item_info">',
#         html)
#     blist = re.findall(r'<em class="item_source">.+</em>', html)
#     for i in range(len(alist)):
#         print('({})'.format(i + 1), end='')
#         print(blist[i][24:-42], end='于')
#         print(blist[i][-15:-5], end='报道：')
#         print(alist[i][23:-33])
#
#
def _Install():
    import os
    libs = [
        "numpy",
        "matplotlib",
        "pillow",
        "sklearn",
        "requests",
        "jieba",
        "pyperclip",
        "wheel",
        "pyautogui",
        "sympy",
        "pyinstaller",
        "Cpython",
        "wheel",
        "Spyder",
        "opencv-python -i",
        "pandas",
        "wordcloud",
        "ipython",
        "xlrd",
        "pygame"]
    # libs=['pyautogui']

    try:
        for lib in libs:
            os.system("pip install " + lib)
            print("Successful")
    except BaseException:
        print("Failed Somehow")
#
#
# def _HollandRadarDraw():
#     import matplotlib.pyplot as plt
#     import matplotlib
#     import numpy as np
#
#     matplotlib.rcParams['font.family'] = 'SimHei'
#     radar_labels = np.array(['研究型(I)', '艺术型(A)', '社会型(S)',
#                              '企业型(E)', '常规型(C)', '现实型(R)'])
#     data = np.array([[0.40, 0.32, 0.35, 0.30, 0.30, 0.88],
#                      [0.85, 0.35, 0.30, 0.40, 0.40, 0.30],
#                      [0.43, 0.89, 0.30, 0.28, 0.22, 0.30],
#                      [0.30, 0.25, 0.48, 0.85, 0.45, 0.40],
#                      [0.20, 0.38, 0.87, 0.45, 0.32, 0.28],
#                      [0.34, 0.31, 0.38, 0.40, 0.92, 0.28]])  # 数据值
#     data_labels = ('艺术家', '实验员', '工程师', '推销员', '社会工作者', '记事员')
#     angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
#     data = np.concatenate((data, [data[0]]))
#     angles = np.concatenate((angles, [angles[0]]))
#     fig = plt.figure(facecolor="white")
#     plt.subplot(111, polar=True)
#     plt.plot(angles, data, 'o-', linewidth=1, alpha=0.2)
#     plt.fill(angles, data, alpha=0.25)
#     plt.thetagrids(angles * 180 / np.pi, radar_labels, frac=1.2)
#     plt.figtext(0.52, 0.95, '霍兰德人格分析', ha='center', size=20)
#     legend = plt.legend(data_labels, loc=(0.94, 0.80), labelspacing=0.1)
#     plt.setp(legend.get_texts(), fontsize='large')
#     plt.grid(True)
#     plt.savefig('holland_radar.jpg')
#     plt.show()
#
#
# def _玫瑰():
#     # RoseDraw.py
#     import turtle as t
#     # 定义一个曲线绘制函数
#
#     def DegreeCurve(n, r, d=1):
#         for i in range(n):
#             t.left(d)
#             t.circle(r, abs(d))
#
#     # 初始位置设定
#     s = 0.2  # size
#     t.setup(450 * 5 * s, 750 * 5 * s)
#     t.pencolor("black")
#     t.fillcolor("red")
#     t.speed(100)
#     t.penup()
#     t.goto(0, 900 * s)
#     t.pendown()
#     # 绘制花朵形状
#     t.begin_fill()
#     t.circle(200 * s, 30)
#     DegreeCurve(60, 50 * s)
#     t.circle(200 * s, 30)
#     DegreeCurve(4, 100 * s)
#     t.circle(200 * s, 50)
#     DegreeCurve(50, 50 * s)
#     t.circle(350 * s, 65)
#     DegreeCurve(40, 70 * s)
#     t.circle(150 * s, 50)
#     DegreeCurve(20, 50 * s, -1)
#     t.circle(400 * s, 60)
#     DegreeCurve(18, 50 * s)
#     t.fd(250 * s)
#     t.right(150)
#     t.circle(-500 * s, 12)
#     t.left(140)
#     t.circle(550 * s, 110)
#     t.left(27)
#     t.circle(650 * s, 100)
#     t.left(130)
#     t.circle(-300 * s, 20)
#     t.right(123)
#     t.circle(220 * s, 57)
#     t.end_fill()
#     # 绘制花枝形状
#     t.left(120)
#     t.fd(280 * s)
#     t.left(115)
#     t.circle(300 * s, 33)
#     t.left(180)
#     t.circle(-300 * s, 33)
#     DegreeCurve(70, 225 * s, -1)
#     t.circle(350 * s, 104)
#     t.left(90)
#     t.circle(200 * s, 105)
#     t.circle(-500 * s, 63)
#     t.penup()
#     t.goto(170 * s, -30 * s)
#     t.pendown()
#     t.left(160)
#     DegreeCurve(20, 2500 * s)
#     DegreeCurve(220, 250 * s, -1)
#     # 绘制一个绿色叶子
#     t.fillcolor('green')
#     t.penup()
#     t.goto(670 * s, -180 * s)
#     t.pendown()
#     t.right(140)
#     t.begin_fill()
#     t.circle(300 * s, 120)
#     t.left(60)
#     t.circle(300 * s, 120)
#     t.end_fill()
#     t.penup()
#     t.goto(180 * s, -550 * s)
#     t.pendown()
#     t.right(85)
#     t.circle(600 * s, 40)
#     # 绘制另一个绿色叶子
#     t.penup()
#     t.goto(-150 * s, -1000 * s)
#     t.pendown()
#     t.begin_fill()
#     t.rt(120)
#     t.circle(300 * s, 115)
#     t.left(75)
#     t.circle(300 * s, 100)
#     t.end_fill()
#     t.penup()
#     t.goto(430 * s, -1070 * s)
#     t.pendown()
#     t.right(30)
#     t.circle(-600 * s, 35)
#     t.done()
#
#
# def 自动化():
#     import pyautogui
#     import time
#     import xlrd
#     import pyperclip
#
#     # 定义鼠标事件
#
#     # pyautogui库其他用法
#     # https://blog.csdn.net/qingfengxd1/article/details/108270159
#
#     def mouseClick(clickTimes, lOrR, img, reTry):
#         if reTry == 1:
#             while True:
#                 location = pyautogui.locateCenterOnScreen(img, confidence=0.9)
#                 if location is not None:
#                     pyautogui.click(
#                         location.x,
#                         location.y,
#                         clicks=clickTimes,
#                         interval=0.2,
#                         duration=0.2,
#                         button=lOrR)
#                     break
#                 print("未找到匹配图片,0.1秒后重试")
#                 time.sleep(0.1)
#         elif reTry == -1:
#             while True:
#                 location = pyautogui.locateCenterOnScreen(img, confidence=0.9)
#                 if location is not None:
#                     pyautogui.click(
#                         location.x,
#                         location.y,
#                         clicks=clickTimes,
#                         interval=0.2,
#                         duration=0.2,
#                         button=lOrR)
#                 time.sleep(0.1)
#         elif reTry > 1:
#             i = 1
#             while i < reTry + 1:
#                 location = pyautogui.locateCenterOnScreen(img, confidence=0.9)
#                 if location is not None:
#                     pyautogui.click(
#                         location.x,
#                         location.y,
#                         clicks=clickTimes,
#                         interval=0.2,
#                         duration=0.2,
#                         button=lOrR)
#                     print("重复")
#                     i += 1
#                 time.sleep(0.1)
#
#     # 数据检查
#     # cmdType.value  1.0 左键单击    2.0 左键双击  3.0 右键单击  4.0 输入  5.0 等待  6.0 滚轮
#     # ctype     空：0
#     #           字符串：1
#     #           数字：2
#     #           日期：3
#     #           布尔：4
#     #           error：5
#     def dataCheck(sheet1):
#         checkCmd = True
#         # 行数检查
#         if sheet1.nrows < 2:
#             print("无数据")
#             checkCmd = False
#         # 每行数据检查
#         i = 1
#         while i < sheet1.nrows:
#             # 第1列 操作类型检查
#             cmdType = sheet1.row(i)[0]
#             if cmdType.ctype != 2 or (cmdType.value != 1.0 and cmdType.value != 2.0 and cmdType.value != 3.0
#                                       and cmdType.value != 4.0 and cmdType.value != 5.0 and cmdType.value != 6.0):
#                 print('第', i + 1, "行,第1列数据有误")
#                 checkCmd = False
#             # 第2列 内容检查
#             cmdValue = sheet1.row(i)[1]
#             # 读图点击类型指令，内容必须为字符串类型
#             if cmdType.value == 1.0 or cmdType.value == 2.0 or cmdType.value == 3.0:
#                 if cmdValue.ctype != 1:
#                     print('第', i + 1, "行,第2列数据误")
#                     checkCmd = False
#             # 输入类型，内容不能为空
#             if cmdType.value == 4.0:
#                 if cmdValue.ctype == 0:
#                     print('第', i + 1, "行,第2列数据误")
#                     checkCmd = False
#             # 等待类型，内容必须为数字
#             if cmdType.value == 5.0:
#                 if cmdValue.ctype != 2:
#                     print('第', i + 1, "行,第2列数据有误")
#                     checkCmd = False
#             # 滚轮事件，内容必须为数字
#             if cmdType.value == 6.0:
#                 if cmdValue.ctype != 2:
#                     print('第', i + 1, "行,第2列数据有误")
#                     checkCmd = False
#             i += 1
#         return checkCmd
#
#     # 任务
#     def mainWork(img):
#         i = 1
#         while i < sheet1.nrows:
#             # 取本行指令的操作类型
#             cmdType = sheet1.row(i)[0]
#             if cmdType.value == 1.0:
#                 # 取图片名称
#                 img = sheet1.row(i)[1].value
#                 reTry = 1
#                 if sheet1.row(i)[2].ctype == 2 and sheet1.row(i)[2].value != 0:
#                     reTry = sheet1.row(i)[2].value
#                 mouseClick(1, "left", img, reTry)
#                 print("单击左键", img)
#             # 2代表双击左键
#             elif cmdType.value == 2.0:
#                 # 取图片名称
#                 img = sheet1.row(i)[1].value
#                 # 取重试次数
#                 reTry = 1
#                 if sheet1.row(i)[2].ctype == 2 and sheet1.row(i)[2].value != 0:
#                     reTry = sheet1.row(i)[2].value
#                 mouseClick(2, "left", img, reTry)
#                 print("双击左键", img)
#             # 3代表右键
#             elif cmdType.value == 3.0:
#                 # 取图片名称
#                 img = sheet1.row(i)[1].value
#                 # 取重试次数
#                 reTry = 1
#                 if sheet1.row(i)[2].ctype == 2 and sheet1.row(i)[2].value != 0:
#                     reTry = sheet1.row(i)[2].value
#                 mouseClick(1, "right", img, reTry)
#                 print("右键", img)
#             # 4代表输入
#             elif cmdType.value == 4.0:
#                 inputValue = sheet1.row(i)[1].value
#                 pyperclip.copy(inputValue)
#                 pyautogui.hotkey('ctrl', 'v')
#                 time.sleep(0.5)
#                 print("输入:", inputValue)
#             # 5代表等待
#             elif cmdType.value == 5.0:
#                 # 取图片名称
#                 waitTime = sheet1.row(i)[1].value
#                 time.sleep(waitTime)
#                 print("等待", waitTime, "秒")
#             # 6代表滚轮
#             elif cmdType.value == 6.0:
#                 # 取图片名称
#                 scroll = sheet1.row(i)[1].value
#                 pyautogui.scroll(int(scroll))
#                 print("滚轮滑动", int(scroll), "距离")
#             i += 1
#
#     if __name__ == '__main__':
#         # file = 'cmd.xls'
#         # 打开文件
#         file = input('输入路径：')
#         wb = xlrd.open_workbook(filename=file)
#         # 通过索引获取表格sheet页
#         sheet1 = wb.sheet_by_index(0)
#         print('开始')
#         # 数据检查
#         checkCmd = dataCheck(sheet1)
#         if checkCmd:
#             key = input('选择功能: 1.做一次 2.死循环 \n')
#             if key == '1':
#                 # 循环拿出每一行指令
#                 mainWork(sheet1)
#             elif key == '2':
#                 while True:
#                     mainWork(sheet1)
#                     time.sleep(0.1)
#                     print("等待0.1秒")
#         else:
#             print('输入有误或者已经退出!')
#
#
# def 自动绘制():
#     import turtle as t
#
#     # D:/Programming/Python/txt/data.txt
#     def draw():
#         t.title("自动绘制")
#         t.setup(800, 600)
#         t.pensize(6)
#         t.pencolor('red')
#         file = open(input())
#         t.hideturtle()
#         t.print(file)
#         data = []
#         for line in file:
#             line = line.replace('\n', '0')
#             data.append(list(map(eval, line.split(','))))
#         file.close()
#         for i in range(len(data) - 1):
#             t.pencolor(data[i][3], data[i][4], data[i][5] % 1)
#             t.fd(data[i][0])
#             if data[i][1]:
#                 t.right(data[i][2])
#             else:
#                 t.left(data[i][2])
#         t.exitonclick()
#
#
# def _IP():
#     import requests
#     import os
#
#     def GetPhotoFromWeb():
#         root = "C:\\Users\\15149\\Desktop\\"
#         url2 = input('the position of file:')
#         vk = {'user-agent': 'Edge/97.0.1072.62'}
#         # url2='http://cms-bucket.ws.126.net/2022 /0113/3db22776j00r5n4y4003cc000ow009vc.jpg'
#         path = root + url2.split('/')[-1]
#         try:
#             if not os.path.exists(root):
#                 os.mkdir(root)
#             if not os.path.exists(path):
#                 r = requests.get(url2, params=vk)
#                 with open(path, 'wb') as f:
#                     f.write(r.content)
#                     f.close()
#                     print('successful')
#             else:
#                 print('the file has been existed')
#         except BaseException:
#             print('fail')
#
#     def GetIpPosition():
#         try:
#             ip = input("input :(手机号加db/;ip加ip/)")
#             url = 'https://ip.cn/'
#             suf = '.html'
#             vk = {'user-agent': 'Edge/97.0.1072.62'}
#             r = requests.get(url + ip + suf, params=vk)
#             r.raise_for_status()
#             r.encoding = r.apparent_encoding
#             print(r.text[3600:4200])
#         except BaseException:
#             print("failure")
#
#     # 202.204.80.112
#
#     # GetPhotoFromWeb()
#     GetIpPosition()
#
#
# def _学校排名():
#     import bs4.element
#     import re
#     import requests
#
#     def GetHtmlText(url):
#         try:
#             r = requests.get(url, timeout=30)
#             r.raise_for_status()
#             r.encoding = r.apparent_encoding
#             return r.text
#         except BaseException:
#             return ''
#
#     def FillUnivList(ulist, html):
#         soup = bs4.BeautifulSoup(html, 'html.parser')
#         for i in soup.find_all(True, re.compile('name-cn')):
#             print(i.string)
#
#     def PrintUnivList(ulist, num):
#         tplt = "{0:^10}\t{1:{3}^10}\t{2:^10}"
#         print(tplt.format('排名', '学校名称', '总分', chr(12288)))
#         for i in range(num):
#             u = ulist[i]
#             print(tplt.format(u[0], u[1], u[2], chr(12288)))
#
#     def main():
#         uinfo = []
#         url = 'https://www.shanghairanking.cn/rankings/bcsr/2021/0827'
#         html = GetHtmlText(url)
#         FillUnivList(uinfo, html)
#         # PrintUnivList(uinfo,10)
#
#     main()
#
#
# def _圆周率1():
#     from random import random
#     from time import perf_counter
#     hit = 0
#     N = int(input('input:'))
#     start = perf_counter()
#     for i in range(1, N + 1):
#         x, y = random(), random()
#         if x * x + y * y <= 1.0:
#             hit += 1
#     print("pi={},time={:.5f}s".format(hit / N * 4, perf_counter() - start))
#
#
# def _淘宝():
#     # CrowTaobaoPrice.py
#     import requests
#     import re
#
#     def getHTMLText(url):
#         try:
#             vk = {'useragent': 'Edge/ 97.0.1072.62'}
#             r = requests.get(url, timeout=30, params=vk)
#             r.raise_for_status()
#             r.encoding = r.apparent_encoding
#             # print(r.headers)
#             print(r.status_code)
#             print(r.text)
#             return r.text
#         except BaseException:
#             return ""
#
#     def parsePage(ilt, html):
#         try:
#             # print(plt)
#             # del plt
#             # plt = re.findall(r'\"view_price\"\:\"[\d\.]*\"', html)
#
#             plt = re.findall(r'"view_price":"\d*\.\d*?"', html)
#             # plt=re.findall(compile('"view_price":"\d*.\d*"'))
#             tlt = re.findall(r'"raw_title":".*?"', html)
#             for i in range(len(plt)):
#                 price = eval(plt[i].split(':')[1])
#                 title = eval(tlt[i].split(':')[1])
#                 ilt.append([price, title])
#         except BaseException:
#             print("")
#
#     def printGoodsList(ilt):
#         tplt = "{:4}\t{:8}\t{:16}"
#         print(tplt.format("序号", "价格", "商品名称"))
#         count = 0
#         for g in ilt:
#             count = count + 1
#             print(tplt.format(count, g[0], g[1]))
#
#     def main():
#         goods = '书包'
#         depth = 3
#         start_url = 'https://s.taobao.com/search?initiative_id=staobaoz_20220118&q=' + goods
#         infoList = []
#         for i in range(depth):
#             try:
#                 url = start_url + \
#                     '&bcoffset=1&p4ppushleft=2%2C48&ntoffset=1&s=' + str(44 * i)
#                 html = getHTMLText(url)
#                 # print(html)
#                 parsePage(infoList, html)
#             except BaseException:
#                 continue
#         printGoodsList(infoList)
#
#     # main()
#     def main2():
#         infoList = []
#         getHTMLText
#         f = open('D:/Programming/Python/txt/taobao.txt', 'r', encoding='utf-8')
#         html = f.read()
#         f.close()
#         parsePage(infoList, html)
#         # print(infoList)
#         printGoodsList(infoList)
#
#     main2()
#
#
# def _圆周率2():
#     def f(x):
#         return (4 / (8 * x + 1) - 2 / (8 * x + 4) - 1 /
#                 (8 * x + 5) - 1 / (8 * x + 6)) * pow(16, -x)
#
#     a = 0
#     for i in range(100):
#         a += f(i)
#     print(a)
#
#
# def _daydayup(rate):
#     day = 1
#     for i in range(365):
#         if i % 7 in [0, 6]:
#             day *= 0.99
#         else:
#             day *= (rate + 1)
#     return day
#
#
# def _英文词频():
#     def get_text():
#         print("先下载到同一文件夹再输入英文文件名(带后缀)")
#         txt = open(input(), "r").read().lower()
#         use = [' ']
#         for i in range(26):
#             use.append(chr(ord('a') + i))
#         for i in txt:
#             if i not in use:
#                 txt = txt.replace(i, ' ')
#         return txt
#
#     # "D:/Programming/Python/dist/hamlet.txt"
#
#     def main():
#         txt = get_text()
#         words = txt.split()
#         count = {}
#         for word in words:
#             count[word] = count.get(word, 0) + 1
#         t = list(count.items())
#         t.sort(key=lambda x: x[1], reverse=True)
#         for i in range(10):
#             word, count = t[i]
#             print("第{}名->{}\n{:<6}次".format(i + 1, word, count))
#         input()
#
#
# def _中文词频():
#     import jieba
#
#     def get_text():
#         print("先下载到同一文件夹再输入中文文件名(带后缀)")
#         txt = open(input(), "r", encoding="utf-8").read()
#         txt = jieba.lcut(txt)
#         return txt
#
#     def main():
#         words = get_text()
#         count = {}
#         for word in words:
#             if len(word) == 1:
#                 continue
#             else:
#                 count[word] = count.get(word, 0) + 1
#         t = list(count.items())
#         t.sort(key=lambda x: x[1], reverse=True)
#         for i in range(10):
#             word, count = t[i]
#             print("第{}名->{}\n{:<6}次".format(i + 1, word, count))
#         input()
#
#     main()
#
#
# def _比赛模拟():
#     from random import random
#
#     def PrintIntro():
#         print('这个程序模拟的是两个选手A和B之间的某种比赛(15轮一场)')
#         print('这个模拟需要两个选手的能力值（0~1）')
#
#     def GetInput():
#         a = eval(input('A选手的能力值:'))
#         b = eval(input('B选手的能力值:'))
#         n = eval(input('比赛的场次:'))
#         return a, b, n
#
#     def PrintSummary(WinA, WinB):
#         n = WinA + WinB
#         # print(WinA)
#         # print(WinB)
#         print("现在A与B的{}场模拟开始".format(n))
#         print('在{}场比赛中A获胜{}场，胜率为{:0.1%}'.format(n, WinA, WinA / n))
#         print('在{}场比赛中B获胜{}场，胜率为{:0.1%}'.format(n, WinB, WinB / n))
#
#     def GameOver(wina, winb):
#         return wina == 15 or winb == 15
#
#     def SimOneGame(a, b):
#         wina, winb = 0, 0
#         if random() > 0.5:
#             first = 'A'
#         else:
#             first = 'B'
#         while not GameOver(wina, winb):
#             if first == 'A':
#                 if random() <= a:
#                     wina += 1
#                 first = 'B'
#             else:
#                 if random() <= b:
#                     winb += 1
#                 first = 'A'
#             # print(wina,winb)
#         return wina > winb
#
#     def SimNGame(a, b, n):
#         WinA, WinB = 0, 0
#         for i in range(n):
#             if SimOneGame(a, b):
#                 WinA += 1
#             else:
#                 WinB += 1
#         return WinA, WinB
#
#     def main():
#         PrintIntro()
#         a, b, n = GetInput()
#         WinA, WinB = SimNGame(a, b, n)
#         PrintSummary(WinA, WinB)
#
#     main()
#
#
# def _手绘风格():
#     from PIL import Image
#     import numpy as np
#
#     a = np.asarray(Image.open(input('input(/):')).convert('L')).astype('float')
#
#     depth = 10.  # (0-100)
#     grad = np.gradient(a)  # 取图像灰度的梯度值
#     grad_x, grad_y = grad  # 分别取横纵图像梯度值
#     grad_x = grad_x * depth / 100.
#     grad_y = grad_y * depth / 100.
#     A = np.sqrt(grad_x ** 2 + grad_y ** 2 + 1.)
#     uni_x = grad_x / A
#     uni_y = grad_y / A
#     uni_z = 1. / A
#
#     vec_el = np.pi / 2.2  # 光源的俯视角度，弧度值
#     vec_az = np.pi / 4.  # 光源的方位角度，弧度值
#     dx = np.cos(vec_el) * np.cos(vec_az)  # 光源对x 轴的影响
#     dy = np.cos(vec_el) * np.sin(vec_az)  # 光源对y 轴的影响
#     dz = np.sin(vec_el)  # 光源对z 轴的影响
#
#     b = 255 * (dx * uni_x + dy * uni_y + dz * uni_z)  # 光源归一化
#     b = b.clip(0, 255)
#
#     im = Image.fromarray(b.astype('uint8'))  # 重构图像
#     im.save('D:/Programming/Python/ico&png&jpg/zx1.png')
#
#
# def grey():
#     import numpy
#     from PIL import Image
#     p = input('input:(/)')
#     a = numpy.asarray(Image.open(p).convert('L'))
#     im = Image.fromarray(a.astype('uint8'))
#     im.save(r'D:\Programming\Python\ico&png&jpg\zx3.png')
#
#
# def _中文词云():
#     import jieba
#     import wordcloud
#
#     def ChineseWordcloud():
#         try:
#             f = open(
#                 input("输入中文文件路径（“/”分隔）(D:/Programming/Python/txt/threekingdoms.txt)\n"),
#                 "r",
#                 encoding="utf-8")
#             t = f.read()
#             f.close()
#             ls = jieba.lcut(t)
#             for i in ls:
#                 if len(i) < 2:
#                     del ls[ls.index(i)]
#             txt = " ".join(ls)
#             w = wordcloud.WordCloud(
#                 font_path="msyh.ttc",
#                 width=2000,
#                 height=1500,
#                 background_color="white")
#             w.generate(txt)
#             w.to_file("D:/Programming/Python/ico&png&jpg/wordcloud.png")
#         except BaseException:
#             print("输入错误")
#             ChineseWordcloud()
#
#     ChineseWordcloud()
#
#
# def _七段数码管(x):
#     x = input('input:')
#     import turtle as _
#     t = 0
#     a = {
#         0: (
#             0, 1, 1, 1, 1, 1, 1, 0), 1: (
#             0, 0, 0, 1, 1, 0, 0, 0), 2: (
#             1, 0, 1, 1, 0, 1, 1, 0), 3: (
#             1, 1, 1, 0, 0, 1, 1, 0), 4: (
#             1, 1, 0, 0, 1, 0, 1, 0), 5: (
#             1, 1, 1, 0, 1, 1, 0, 0), 6: (
#             1, 1, 1, 1, 1, 1, 0, 0), 7: (
#             0, 1, 0, 0, 0, 1, 1, 0), 8: (
#             1, 1, 1, 1, 1, 1, 1, 0), 9: (
#             1, 1, 1, 0, 1, 1, 1, 0), }
#     _.setup(1920, 1080)
#     _.color('blue')
#     _.penup()
#     _.pensize(10)
#     _.hideturtle()
#     _.speed(100)
#     _.fd(-940)
#     for i in x:
#         for j in a[int(i)]:
#             _.pendown() if j else _.penup()
#             _.fd(150)
#             _.right(90)
#             t += 1
#             if t % 4 == 0:
#                 _.left(90)
#             if t % 8 == 0:
#                 _.right(180)
#                 _.fd(170)
#     _.exitonclick()
#
#
# def _雪花描绘():
#     import turtle as _
#
#     def _6_3(size, n):
#         if n == 0:
#             _.fd(size)
#         else:
#             for angle in [0, 60, -120, 60]:
#                 _.left(angle)
#                 _6_3(size / 3, n - 1)
#
#     def _6_4():
#         global temporary_size
#         global temporary_n
#         try:
#             temporary_size = int(input())
#             temporary_n = int(input())
#         except ValueError:
#             print("输入有误,请再输一遍")
#             _6_4()
#
#     def main():
#         _6_4()
#         size, n = temporary_size, temporary_n
#         _.setup(1920, 1080)
#         _.hideturtle()
#         _.penup()
#         _.goto(-size / 2, size / (2 * (3 ** 0.5)))
#         _.pendown()
#         _.speed(0)
#         _.pensize(2)
#         _.pencolor('blue')
#         _6_3(size, n)
#         _.right(120)
#         _6_3(size, n)
#         _.right(120)
#         _6_3(size, n)
#         print("点击雪花退出")
#         _.exitonclick()
#
#     print("输入大小与阶数,大小最好小于1000,enter输入")
#     main()
#
#
# def _进度条():
#     import time
#     scale = 50
#     print("执行开始".center(scale + 14, "_"))
#     start = time.perf_counter()
#     for i in range(scale + 1):
#         print("\r{:^3.0f}%[{}->{}]{:.2f}s".format((i / scale) * 100, '*' *
#               i, '..' * (scale - i), time.perf_counter() - start), end=" ")
#         time.sleep(0.05)
#     print("\n" + "执行结束".center(scale + 14, "_"))
#
#
# def _画蛇():
#     import turtle as _
#
#     _.setup(1000, 350, 400, 200)
#     _.penup()
#     _.fd(-250)
#     _.pendown()
#     _.pensize(25)
#     _.pencolor("yellow")
#     _.seth(-40)
#     for i in range(4):
#         _.circle(40, 80)
#         _.circle(-40, 80)
#     _.circle(40, 80 / 2)
#     _.fd(40)
#     _.circle(16, 180)
#     _.fd(40 * 2 / 3)
#     _.done()
#
#
# def _fib0(n):
#     a, b = 0, 1
#     for i in range(n):
#         print("\t", b)
#         a, b = b, a + b
#
#
# def _fib1(n):
#     a, b = 0, 1
#     while b < n:
#         print(b, end=' ')
#         a, b = b, a + b
#     print()
#
#
# def _fib2(n):
#     result = []
#     a, b = 0, 1
#     while b < n:
#         result.append(b)
#         a, b = b, a + b
#     print(result)
#     return result
#
#
# def _hanoi(n=10, a='a', c='c', b='b'):
#     global count
#     if n == 1:
#         print("{}:{}->{}".format(1, a, c))
#         count += 1
#         print(count)
#     else:
#         _hanoi(n - 1, a, b, c)
#         print("{}:{}->{}".format(n, a, c))
#         count += 1
#         print(count)
#         _hanoi(n - 1, b, c, a)
