import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.integrate as integrate
import mpmath as mp
import Exact_Ising_model
import os
import csv

def csv_statistic(data):
    return data.mean(axis='columns'), data.std(axis='columns')
