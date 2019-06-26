import zipfile




with open('./GA7411_TGACCA_L007_R1_002_fastqc.zip') as read:
    content=read.readlines
    print content[0]
    print content[-1]