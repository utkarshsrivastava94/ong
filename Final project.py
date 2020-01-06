#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn; seaborn.set()
import datetime as dt
from sklearn.metrics import r2_score as r2
from datetime import date
from datetime import datetime
from pandas import DataFrame as df
from statsmodels.tsa.arima_model import ARIMA
from pmdarima.arima import auto_arima
from statsmodels.tsa.stattools import adfuller
get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


welldata = pd.read_csv('Daily_Production.csv')
welldata['Date'] = pd.to_datetime(welldata['Date'], format = "%m/%d/%Y %H:%M:%S")
welldata.dtypes
welldata['Date'] = pd.to_datetime(welldata['Date'], format='%m/%d/%Y')
prod_data = welldata[['API', 'Date', 'Production', 'Line_Press']].copy()
prod_data = prod_data.dropna()
prod_data = prod_data[prod_data['Production'] != 0]


# In[ ]:


welldata.dtypes 
prod_data.dtypes
# prod_data.describe()
# welldata['Date']


# In[ ]:


w_49035224030000 = prod_data[prod_data['API'] == 49035224030000]
w_49035224030000 = w_49035224030000[w_49035224030000['Production'] != 0]
plt.plot(w_49035224030000['Date'], w_49035224030000['Production'])
plt.xlabel('Date (years)')
plt.ylabel('Gas Production (Mscf/day)')
# plt.title('Production curve of gas producing well')
plt.show()


# In[ ]:


#prod_data.dtypes
per_well = prod_data.groupby(['API'])['Production'].sum()
#per_well.columns = ['API', 'Production']
count = per_well.min()

p50 = (np.percentile(per_well, 50))
p50 = np.argwhere(per_well == p50).ravel()
well50 = per_well.iloc[p50].index[0]

p10 = (np.percentile(per_well, 10))
p10 = abs(per_well - p10).values.argmin().ravel()
well10 = per_well.iloc[p10].index[0]

p90 = (np.percentile(per_well, 90))
p90 = abs(per_well - p90).values.argmin().ravel()
well90 = per_well.iloc[p90].index[0]

well = []
well.append(well10)
well.append(well50)
well.append(well90)


# In[ ]:


# perc = per_well.iloc[p50].index[0]
# perc
count
per_well.plot()
plt.ylabel('Total production in Mscf per well')
plt.show()
count


# In[ ]:


g = prod_data.groupby(['API'])


well_name = {}
for key in g.indices.keys():
    well_name[key] = g.get_group(key)
    well_name[key] = well_name[key][well_name[key]['Production'] != 0]
    well_name[key]['logProduction'] = np.log(well_name[key]['Production'])
    well_name[key].set_index('Date')
    well_name[key]['Date'] = pd.to_datetime(well_name[key]['Date'], format='%m/%d/%Y')
    
result = adfuller(well_name[well[0]]['logProduction'])
print('ADF Statistic: %f' % result[0])
print('p-value: %f' % result[1])
print('Critical Values:')
for key, value in result[4].items():
	print('\t%s: %.3f' % (key, value))


# In[ ]:


# g.head()
# well_name[49035217710000].set_index('Date')
# well_name[49035217710000]['Date'] = pd.to_datetime(well_name[49035217710000]['Date'])


# In[ ]:


# plt.plot(well_name[49035217710000]['Date'],well_name[49035217710000]['Production'])
for i in well:
    well_name[i].plot(x = 'Date', y ='Production')
    plt.xlabel('Date (years)')
    plt.ylabel('Gas Production (Mscf/day)')
    plt.title('Production curve of gas producing well (API = %d)' % i)
    plt.show()


# In[ ]:


# well_name[49035217710000]['logProduction']= well_name[49035217710000]['logProduction'].astype('float')
# model = ARIMA(well_name[49035217710000]['logProduction'], order=(3,1,3))
# model_fit = model.fit(disp=1)
# print(model_fit.summary())


# In[ ]:


# X = well_name[49035217710000]['Production']
# size = int(len(X) * 0.66)
# train, test = X[0:size], X[size:len(X)]
# history = [x for x in train]
# predictions = list()
# for t in range(len(test)):
# 	model = ARIMA(history, order=(3,1,3))
# 	model_fit = model.fit(disp=0)
# 	output = model_fit.forecast()
# 	yhat = output[0]
# 	predictions.append(yhat)


# In[ ]:


# model_fit.forecast()


# In[ ]:


# fig, ax = plt.subplots(1,1,figsize=(15,5))
# index=np.linspace(size+1,size+1+np.shape(test)[0],np.shape(test)[0])
# plt.plot(index,test,'.g',label='test data')
# plt.plot(index,predictions, 'or', label='predictions')
# plt.plot(train, 'k-', label='decline curve')
# plt.legend(loc=0)
# plt.xlabel('Month')
# plt.ylabel('Gas Production (Mscf/day)')
# plt.title('Production graph showing predicted data using ARIMA model')
# plt.show()


# In[ ]:


import pmdarima as pm

model = pm.auto_arima(well_name[well[2]]['Production'], start_p=1, start_q=1,
#                       test = 'adf',    # using adftest to find optimal 'd'
                      max_p=3, max_q=3, # maximum p and q
                      m=1,              # frequency of series
                      d = 1,           # letting model determine 'd'
                      seasonal=False,   # Assuming Seasonality
                      start_P=0, 
                      trace=True,
                      D=0,
                      error_action='ignore',  
                      suppress_warnings=True, 
                      stepwise=True)
    
print('Model for well with API = %d' % well[2])
print(model.summary())


# In[ ]:


model.plot_diagnostics(figsize=(15,15))
plt.show()


# In[ ]:


# Forecast
n_periods = 180
fc, confint = model.predict(n_periods=n_periods, return_conf_int=True)
# index_of_fc = np.arange(len(well_name[49035217710000]['logProduction']), len(well_name[49035217710000]['logProduction'])+n_periods)
index_of_fc = pd.date_range(start = well_name[well[2]]['Date'].iloc[-n_periods], periods=n_periods)
# make series for plotting purpose
fc_series = pd.Series(fc, index=index_of_fc)
#fc_series.index = pd.to_datetime(fc_series.index)
#lower_series = pd.Series(confint[:, 0], index=index_of_fc)
#upper_series = pd.Series(confint[:, 1], index=index_of_fc)

# Plot
well_name[well[2]].plot(x = 'Date', y ='Production', label = 'Actual production')
fc_series.plot(label='Predicted production')
# plt.fill_between(lower_series.index, 
#                    lower_series, 
#                    upper_series, 
#                    color='k', alpha=.15)
plt.ylabel('Gas production (Mscf)')
plt.title("Final Forecast of Production")
plt.legend()
plt.show()


# In[ ]:


well_name[well[2]]['Production'].iloc[-n_periods:]

# def rsquared(x, y):
#     """ Return R^2 where x and y are array-like."""

#     slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
#     return r_value**2

# rsquared(fc, well_name[well[2]].iloc[-n_periods:])


# In[ ]:


RMSE = ((fc - well_name[well[2]]['Production'].iloc[-n_periods:]) ** 2).mean() ** .5
print('RMSE for P90 well is %d' % RMSE)


# In[ ]:


ME = (fc.mean() - well_name[well[2]]['Production'].iloc[-n_periods:].mean())
ME

