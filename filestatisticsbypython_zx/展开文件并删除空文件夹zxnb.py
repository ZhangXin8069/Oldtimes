def fileStatisticsbyPython_zx(pah):
    from os import listdir
    from os.path import isdir, join, split
    from pandas import DataFrame
    lit = []
    for i in listdir(pah):
        if isdir(join(pah, i)):
            lit.append(join(pah, i))
    print(lit)
    lc = []
    for i in lit:
        lt = list('=HYPERLINK(\"' + join(pah, i, j) +
                  '\",\"' + j + '\")' for j in listdir(join(pah, i)))
        lt.insert(0, str(len(listdir(join(pah, i)))))
        lc.append(lt)
    df = DataFrame(data=lc, index=lit)
    df.to_excel(join(pah, 'fileStatisticsbyPython_zx.xlsx'), na_rep=' ')
    print(df)


def findAllfiles(pah):
    from os import walk, path
    lit = []
    for foldName, subfolders, filenames in walk(pah):
        for filename in filenames:
            try:
                lit.append(path.join(foldName, filename))
            except Exception as e:
                #print(e)
                pass
    print(lit)
    return lit


def fileSort(pah=[], want_lit=[], file_lit=[]):
    from os import mkdir, rename
    from os.path import join, split
    for i in want_lit:
        try:
            mkdir(join(pah, i))
        except Exception as e:
            #print(e)
            pass
        finally:
            for j in file_lit:
                (_, j_) = split(j)
                if i in j_:
                    rename(j, join(pah, i, j_))
def fileClean(pah):
    from os import walk,rmdir
    from os.path import join
    for (root, dirs, files) in walk(pah):
        for item in dirs:
            dir = join(root, item)
            try:
                rmdir(dir)  #os.rmdir() 方法用于删除指定路径的目录。仅当这文件夹是空的才可以, 否则, 抛出OSError。
                print(dir)
            except Exception as e:
                print('Exception',e)

if __name__ == '__main__':
    want_lit=['郭宇飞','贺睿康','贺行伟','黄杰峰','黄啸','蒋永亮','康博生','马钰卓','石芙蓉','孙政','王国栋','易佳','张鑫','张阳洁','钟杰']
    want_lit=['']
    pah = 'E:/desktop/21核物二班核酸检测记录截图'
    # absolute path
    lit = findAllfiles(pah=pah)
    fileSort(file_lit=lit, pah=pah,want_lit=want_lit)
    fileStatisticsbyPython_zx(pah=pah)
    fileClean(pah=pah)