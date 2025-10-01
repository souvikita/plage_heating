import museval
from museval.utils import save_eis_iris_dates
import os
from tqdm import tqdm

urls = ['https://www.lmsal.com/hek/hcr?cmd=search-events-corr&outputformat=json&startTime=2017-07-29T23:00&stopTime=2017-07-30T02:00&instrument=IRIS&hasData=true&hideMostLimbScans=true&optionalcorr=SOT&optionalcorr=SOTSP&optionalcorr=XRT&requiredcorr=EIS']
text_files_path = '/Users/souvikb/various_analysis/GaussSep/bose_codes/'
for ii, url in tqdm(enumerate(urls)):
    save_eis_iris_dates([url], f"{text_files_path}HOP_307_dates_{ii}.txt", alternate_only=False)
    