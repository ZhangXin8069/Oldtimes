def news():
    import requests
    import re
    print('国际新闻')
    url = 'https://news.china.com/international/'
    r = requests.get(url)
    r.encoding = 'utf-8'
    html = r.text
    alist = re.findall(
        r'.html" target="_blank">.+</a></h3><span class="item_info">',
        html)
    blist = re.findall(r'<em class="item_source">.+</em>', html)
    for i in range(20):
        print('({})'.format(i + 1), end='')
        print(blist[i][24:-42], end='于')
        print(blist[i][-15:-5], end='报道：')
        print(alist[i][23:-33])
    print('兰大新闻')
    url = 'http://jwc.lzu.edu.cn/'
    r = requests.get(url)
    r.encoding = 'GB2312'
    html = r.text
    li = re.findall(
        r'<li><a title=\'(.+?)\'  href=\".*?/lzupage(.+?)\"',
        html)
    for i in range(-10, 0):
        print('({})'.format(i + 11), end='')
        print(li[i][0] + '...http://jwc.lzu.edu.cn/lzupage' + li[i][1])
    input('input for leave')


if __name__ == '__main__':
    news()
