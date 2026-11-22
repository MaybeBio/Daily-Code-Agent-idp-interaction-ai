# Convert amino acid sequence to 3Di sequence using Foldseek
import requests
import time


def _split_seq(seq, window=1000):
    return [seq[(i*window):min((i+1)*window,len(seq))] for i in range((len(seq)-1)//window+1)]


def get_saprotseq(sequence: str):
    if len(sequence) <= 1200:
        # if sequence is short, get saprotseq for the sequence by foldseek api
        response = requests.get(f"https://3di.foldseek.com/predict/{sequence}")
        s2 = response.json()
        time.sleep(0.5)
        saprotseq = ''.join([a+b for a, b in zip(sequence, s2.lower())])
    elif len(sequence) > 1200:
        # if sequence is long, break it down
        # get saprotseq for the subsequences by foldseek api
        s2 = ""
        for ss in _split_seq(sequence):
            response = requests.get(f"https://3di.foldseek.com/predict/{ss}")
            s2 += response.json()
            time.sleep(0.5)
        saprotseq = ''.join([a+b for a, b in zip(sequence, s2.lower())])
    return saprotseq
