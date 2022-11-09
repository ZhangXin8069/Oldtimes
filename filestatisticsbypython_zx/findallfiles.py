def findAllfiles(pa="D:\\Monitor"):
    from os import walk, path, getcwd
    pa = input('path:\n')
    if pa == '':
        pa = getcwd()
    f = open('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt', 'w')
    for foldName, subfolders, filenames in walk(pa):
        for filename in filenames:
            try:
                f.write(path.join(foldName, filename) + '\r')
            except Exception as e:
                print(e)
    f.close()
    print('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt')


def findallfiles(pa="D:\\Monitor"):
    from os import walk, path
    li = []
    f = open('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt', 'w')
    for foldName, subfolders, filenames in walk(pa):
        for filename in filenames:
            try:
                li.append(path.join(foldName, filename))
            except Exception as e:
                print(e)
                pass

    f.close()
    return li
    # print('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt')


if __name__ == '__main__':
    # from os import getcwd
    #
    # pa = r"D:\Programming"
    # print(getcwd())
    # f = open(getcwd() + '\\files.txt', 'w+')
    # for i in findAllfiles(pa):
    #     print(i)
    #     f.write(i + '\r')
    # f.close()
    print(findallfiles())
    # pa = "D:\\"
    # from os import path
    # print(path.getatime(pa))
    # print(path.getmtime(pa))
