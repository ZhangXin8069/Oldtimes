def fileStatisticsbyPython_zx(pa=r'D:\Monitor'):
    from os import listdir, path
    from pandas import DataFrame
    li = []
    for i in listdir(pa):
        if path.isdir(pa + '\\' + i):
            li.append(i)
    lc = []
    for i in li:
        lt = list('=HYPERLINK(\"' + pa + '\\' + i + '\\' + j + '\",\"' + j + '\")' for j in
                  listdir(pa + '\\' + i))
        for k in lt:
            print(k)
            if k[-5:-2] == 'cfg':
                print('ok')
                lt.remove(k)
        lt.insert(0, str(len(lt)))
        lc.append(lt)
    df = DataFrame(data=lc, index=li)
    df.to_excel(pa + '\\' + 'fileStatisticsbyPython_zx.xlsx', na_rep=' ')
    print(df)
    input('input for leave')

from sys import argv
fileStatisticsbyPython_zx(argv[1])
