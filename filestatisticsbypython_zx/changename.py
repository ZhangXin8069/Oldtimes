# def copy(pa='D:\\Monitor', spa='D:\\test'):
#     from shutil import copytree, rmtree
#     rmtree(spa)
#     copytree(pa, spa)
#     from findallfiles import findallfiles
#     print(findallfiles(pa=spa))


def changename(pa='D:\\test2', spa='D:\\test'):
    from os import walk, rename, mkdir
    from os.path import join
    try:
        mkdir(spa)
    except FileExistsError:
        pass
    for foldName, subfolders, filenames in walk(pa):
        # print('foldName:{}'.format(foldName))#全部文件夹中的一个
        # print('subfolders:{}'.format(subfolders))#本文件夹的子文件夹
        # print('filenames::{}'.format(filenames))#本文件夹中的文件
        for filename in filenames:
            abp = join(foldName, filename)
            nen = join(spa, filename)
            try:
                rename(abp, nen)
            except FileExistsError:
                # print(Exception)
                pass


if __name__ == '__main__':
    # copy()
    changename()
