def shu1():
    import shutil
    # from os import mkdir
    #
    # mkdir('D:\\Programming\\learning')
    # shutil.copytree("D:\\Programming\\User\\learning", "D:\\Programming\\learning")
    # shutil.copytree()
    # try:
    #     shutil.move('D:\\Programming\\learning', r"C:\Users\15149\Desktop\新建文件夹")
    # except Exception as e:
    #     print(e)
    # shutil.rmtree("D:\\Programming\\learning")
    # os.unlink(r"D:\Programming\User\learning\.._.txt")
    # shutil.make_archive(r'压缩包', 'zip', root_dir=r"D:\Programming\User\learning\文件处理")
    shutil.unpack_archive(r"D:\Programming\User\filestatisticsbypython_zx\压缩包.zip", format='zip')


def shu2():
    from findallfiles import findAllfiles
    from os import stat
    pa = "D:\\Monitor"
    findAllfiles(pa)

    f = open('_'.join(pa.replace(':', '').split(sep='\\')) + '.txt', 'r')
    for i in f.readlines():
        print(stat(i[:-1]))
    f.close()


def shu3():
    import zipfile
    # z = zipfile.ZipFile(r"D:\Programming\User\filestatisticsbypython_zx\压缩包2.zip", 'w')
    z = zipfile.ZipFile(r"D:\Programming\User\filestatisticsbypython_zx\压缩包2.zip")
    # z.write(r"D:\Programming\User\filestatisticsbypython_zx\D_Programming_User.txt")
    z.extractall(r'D:\Programming\User\filestatisticsbypython_zx')
    z.close()
    # z = zipfile.ZipFile(r"D:\Programming\User\filestatisticsbypython_zx\压缩包2.zip", 'a')
    # z.write(r"D:\Programming\User\zx.ico")
    # z.close()


if __name__ == '__main__':
    shu3()
