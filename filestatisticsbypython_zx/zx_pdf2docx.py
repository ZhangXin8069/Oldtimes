def pdf_to_word(fileName):
    import re

    from pdf2docx import Converter

    pdf_file = fileName
    name = re.findall(r'(.*?)\.', pdf_file)[0]
    docx_file = f'{name}.docx'

    cv = Converter(pdf_file)
    cv.convert(docx_file, start=0, end=None)
    cv.close()
