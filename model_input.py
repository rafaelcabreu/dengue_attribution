import geopandas as gpd
import numpy as np
import pandas as pd
import json

from unidecode import unidecode


# Regions
region_map = json.loads(open('utils/regions.json', 'r').read())


# Temperature
era5_brazil = pd.read_csv('data/era5-tmean-city-brazil-2000-2024.csv', parse_dates=['time'])
era5_brazil['t2m'][era5_brazil['t2m'] < 0] = era5_brazil['t2m'][era5_brazil['t2m'] < 0] + 273.15
era5_brazil = era5_brazil.rename(columns={'t2m': 'tmed'})

era5_mexico = pd.read_csv('data/era5-tmean-city-mexico-2000-2024.csv', parse_dates=['time'])
era5_mexico['tmed'] = era5_mexico['tmed'] - 273.15

era5_peru = pd.read_csv('data/era5-tmean-city-peru-2000-2024.csv', parse_dates=['time'])
era5_peru['tmed'] = era5_peru['tmed'] - 273.15

era5 = pd.concat([era5_brazil[['city_residency', 'time', 'tmed']], era5_mexico[['city_residency', 'time', 'tmed']], era5_peru[['city_residency', 'time', 'tmed']]])
era5 = era5.rename(columns={'time': 'date_first_symptoms'})

tmean0 = era5.pivot_table(index='date_first_symptoms', columns='city_residency', values='tmed').reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='tmed').rename(columns={'tmed': 'mean_2m_air_temp_degree1'})
tmean1 = era5.pivot_table(index='date_first_symptoms', columns='city_residency', values='tmed').shift(1).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='tmed').rename(columns={'tmed': 'mean_2m_air_temp_degree1_lag1'})
tmean2 = era5.pivot_table(index='date_first_symptoms', columns='city_residency', values='tmed').shift(2).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='tmed').rename(columns={'tmed': 'mean_2m_air_temp_degree1_lag2'})
tmean3 = era5.pivot_table(index='date_first_symptoms', columns='city_residency', values='tmed').shift(3).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='tmed').rename(columns={'tmed': 'mean_2m_air_temp_degree1_lag3'})
tmean4 = era5.pivot_table(index='date_first_symptoms', columns='city_residency', values='tmed').shift(4).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='tmed').rename(columns={'tmed': 'mean_2m_air_temp_degree1_lag4'})
tmean5 = era5.pivot_table(index='date_first_symptoms', columns='city_residency', values='tmed').shift(5).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='tmed').rename(columns={'tmed': 'mean_2m_air_temp_degree1_lag5'})

temp = pd.concat([
    tmean0.set_index(['date_first_symptoms', 'city_residency']),
    tmean1.set_index(['date_first_symptoms', 'city_residency']),
    tmean2.set_index(['date_first_symptoms', 'city_residency']),
    tmean3.set_index(['date_first_symptoms', 'city_residency']),
    tmean4.set_index(['date_first_symptoms', 'city_residency']),
    tmean5.set_index(['date_first_symptoms', 'city_residency'])
], axis=1).reset_index()

temp['date_first_symptoms'] = temp['date_first_symptoms'].apply(lambda x: pd.to_datetime(x.strftime('%Y-%m-01')))

temp['mean_2m_air_temp_degree2_lag1'] = temp['mean_2m_air_temp_degree1_lag1'] ** 2
temp['mean_2m_air_temp_degree2_lag2'] = temp['mean_2m_air_temp_degree1_lag2'] ** 2
temp['mean_2m_air_temp_degree2_lag3'] = temp['mean_2m_air_temp_degree1_lag3'] ** 2
temp['mean_2m_air_temp_degree2_lag4'] = temp['mean_2m_air_temp_degree1_lag4'] ** 2
temp['mean_2m_air_temp_degree2_lag5'] = temp['mean_2m_air_temp_degree1_lag5'] ** 2

temp['mean_2m_air_temp_degree3_lag1'] = temp['mean_2m_air_temp_degree1_lag1'] ** 3
temp['mean_2m_air_temp_degree3_lag2'] = temp['mean_2m_air_temp_degree1_lag2'] ** 3
temp['mean_2m_air_temp_degree3_lag3'] = temp['mean_2m_air_temp_degree1_lag3'] ** 3
temp['mean_2m_air_temp_degree3_lag4'] = temp['mean_2m_air_temp_degree1_lag4'] ** 3
temp['mean_2m_air_temp_degree3_lag5'] = temp['mean_2m_air_temp_degree1_lag5'] ** 3

# Precipitation
chirps_brazil = pd.read_csv('data/chirps-precip-monthly-brazil.csv', parse_dates=['time'])
chirps_brazil['state_residency'] = chirps_brazil['city_residency'].apply(lambda x: x.split('_')[0])
chirps_brazil = chirps_brazil[chirps_brazil['time'] >= '2000-01-01']
chirps_brazil['date_first_symptoms'] = chirps_brazil['time'].apply(lambda x: pd.to_datetime(x.strftime('%Y-%m-01')))

chirps_mexico = pd.read_csv('data/chirps-precip-monthly-mexico.csv', parse_dates=['time'])
chirps_mexico = chirps_mexico.groupby('city_residency')[['time', 'precip']].resample('M', on='time').sum().reset_index()
chirps_mexico['date_first_symptoms'] = chirps_mexico['time'].apply(lambda x: pd.to_datetime(x.strftime('%Y-%m-01')))

chirps_peru = pd.read_csv('data/chirps-precip-monthly-peru.csv', parse_dates=['time'])
chirps_peru = chirps_peru.groupby('city_residency')[['time', 'precip']].resample('M', on='time').sum().reset_index()
chirps_peru['date_first_symptoms'] = chirps_peru['time'].apply(lambda x: pd.to_datetime(x.strftime('%Y-%m-01')))

chirps = pd.concat([
    chirps_brazil[['city_residency', 'precip', 'date_first_symptoms']],
    chirps_mexico[['city_residency', 'precip', 'date_first_symptoms']],
    chirps_peru[['city_residency', 'precip', 'date_first_symptoms']]
])

prec0 = chirps.pivot_table(index='date_first_symptoms', columns='city_residency', values='precip').reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='precip').rename(columns={'precip': 'total_precipitation'})
prec1 = chirps.pivot_table(index='date_first_symptoms', columns='city_residency', values='precip').shift(1).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='precip').rename(columns={'precip': 'total_precipitation_lag1'})
prec2 = chirps.pivot_table(index='date_first_symptoms', columns='city_residency', values='precip').shift(2).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='precip').rename(columns={'precip': 'total_precipitation_lag2'})
prec3 = chirps.pivot_table(index='date_first_symptoms', columns='city_residency', values='precip').shift(3).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='precip').rename(columns={'precip': 'total_precipitation_lag3'})
prec4 = chirps.pivot_table(index='date_first_symptoms', columns='city_residency', values='precip').shift(4).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='precip').rename(columns={'precip': 'total_precipitation_lag4'})
prec5 = chirps.pivot_table(index='date_first_symptoms', columns='city_residency', values='precip').shift(5).reset_index().melt(id_vars='date_first_symptoms', var_name='city_residency', value_name='precip').rename(columns={'precip': 'total_precipitation_lag5'})

prec = pd.concat([
    prec0.set_index(['date_first_symptoms', 'city_residency']),
    prec1.set_index(['date_first_symptoms', 'city_residency']),
    prec2.set_index(['date_first_symptoms', 'city_residency']),
    prec3.set_index(['date_first_symptoms', 'city_residency']),
    prec4.set_index(['date_first_symptoms', 'city_residency']),
    prec5.set_index(['date_first_symptoms', 'city_residency'])
], axis=1).reset_index()

# Cases (replicate population from 2020 up to 2024 because DATASUS does not provide it)
cases_brazil = pd.read_csv('data/dengue_city-monthly_brazil.csv', parse_dates=['date_first_symptoms'])
pop = pd.read_csv('data/population_city-yearly_2000-2021.csv')

cases_brazil['year'] = cases_brazil['date_first_symptoms'].dt.year
pop_2021 = pop[pop['year'] == 2020]
pop_2021['year'] = 2021
pop_2022 = pop[pop['year'] == 2020]
pop_2022['year'] = 2022
pop_2023 = pop[pop['year'] == 2020]
pop_2023['year'] = 2023
pop_2024 = pop[pop['year'] == 2020]
pop_2024['year'] = 2024

pop = pd.concat([pop, pop_2021, pop_2022, pop_2023, pop_2024], axis=0)
cases_brazil = cases_brazil.set_index(['city_residency', 'year']).join(pop.set_index(['city_residency', 'year'])['population']).reset_index()
cases_brazil['dengue_inc'] = 1e5 * (cases_brazil['n_cases'] / cases_brazil['population'])
cases_brazil = cases_brazil.dropna().drop_duplicates(subset=['state_residency', 'city_residency', 'date_first_symptoms'])
cases_brazil['city_residency'] = cases_brazil['state_residency'] + '_' + cases_brazil['city_residency']
cases_brazil['country'] = 'Brazil'

cases_mexico = pd.read_csv('data/dengue_city-monthly_mexico.csv', parse_dates=['date_first_symptoms'])
cases_mexico['date_first_symptoms'] = cases_mexico['date_first_symptoms'].apply(lambda x: pd.to_datetime(x.strftime('%Y-%m-01')))
cases_mexico['region'] = cases_mexico['municipio'].map(region_map)
cases_mexico['dengue_inc'] = 1e5 * (cases_mexico['n_cases'] / cases_mexico['population'])
cases_mexico = cases_mexico.rename(columns={'municipio': 'city_residency'})
cases_mexico['country'] = 'Mexico'

cases_peru = pd.read_csv('data/dengue_city-monthly_peru.csv', parse_dates=['date'])
cases_peru['month'] = cases_peru['date'].dt.month
cases_peru['year'] = cases_peru['date'].dt.year
cases_peru = cases_peru.groupby(['region', 'city_residency']).resample('ME', on='date').agg({'n_cases': 'sum', 'POBTOTAL': 'mean'}).reset_index()
cases_peru = cases_peru.rename(columns={'date': 'date_first_symptoms', 'POBTOTAL': 'population'})
cases_peru['date_first_symptoms'] = cases_peru['date_first_symptoms'].apply(lambda x: pd.to_datetime(x.strftime('%Y-%m-01')))
cases_peru['dengue_inc'] = 1e5 * (cases_peru['n_cases'] / cases_peru['population'])
cases_peru['dengue_inc'] = cases_peru['dengue_inc'].fillna(0)
cases_peru['country'] = 'Peru'

cases = pd.concat([
    cases_brazil[['city_residency', 'date_first_symptoms', 'dengue_inc', 'n_cases', 'population', 'country']],
    cases_mexico[['city_residency', 'date_first_symptoms', 'dengue_inc', 'n_cases', 'population', 'country']],
    cases_peru[['city_residency', 'date_first_symptoms', 'dengue_inc', 'n_cases', 'population', 'country']]
])

# mmunity per city
immunity = cases.copy()
immunity['state_residency'] = immunity['city_residency'].apply(lambda x: x.split('_')[0])
immunity['region'] = immunity['state_residency'].map(region_map)
immunity['region'][immunity['country'] == 'Mexico'] = immunity[immunity['country'] == 'Mexico']['city_residency'].map(region_map)

immunity = immunity.groupby(['date_first_symptoms', 'city_residency']).agg({'n_cases': 'sum', 'population': 'sum'}).reset_index()
immunity['year'] = immunity['date_first_symptoms'].apply(lambda x: x.year + 1 if x.month >= 8 else x.year)
immunity = immunity[(immunity['year'] > 2001) & (immunity['year'] < 2025)]
immunity = immunity.groupby(['year', 'city_residency']).agg({'n_cases': 'sum', 'population': 'mean'}).reset_index()

immunity['dengue_inc'] = 1e5 * (immunity['n_cases'] / immunity['population'])
immunity_lag1 = immunity.pivot_table(index='year', columns='city_residency', values='dengue_inc')
immunity_lag1.loc[2025] = np.nan
immunity_lag1 = immunity_lag1.shift(1).reset_index().melt(id_vars='year', var_name='city_residency', value_name='immunity_lag1')

immunity_lag2 = immunity.pivot_table(index='year', columns='city_residency', values='dengue_inc')
immunity_lag2.loc[2025] = np.nan
immunity_lag2.loc[2026] = np.nan
immunity_lag2 = immunity_lag2.shift(2).reset_index().melt(id_vars='year', var_name='city_residency', value_name='immunity_lag2')

immunity_lag3 = immunity.pivot_table(index='year', columns='city_residency', values='dengue_inc')
immunity_lag3.loc[2025] = np.nan
immunity_lag3.loc[2026] = np.nan
immunity_lag3.loc[2027] = np.nan
immunity_lag3 = immunity_lag3.shift(3).reset_index().melt(id_vars='year', var_name='city_residency', value_name='immunity_lag3')

immunity = pd.concat([
    immunity_lag1.set_index(['year', 'city_residency']),
    immunity_lag2.set_index(['year', 'city_residency']),
    immunity_lag3.set_index(['year', 'city_residency']),
], axis=1).reset_index()

# socioeconomic data
ips_data = pd.read_csv('data/ips_brasil_municipios.csv')

ips_data['Município'] = ips_data['Município'].apply(lambda x: str(unidecode(x)).upper()[:-5])
ips_data['city_residency'] = ips_data['UF'] + '_' + ips_data['Município']

columns = ['Água e Saneamento', 'PIB per capita 2021',]
ips_data = ips_data[['city_residency'] + columns].rename(columns={'Água e Saneamento': 'water_sanitation', 'PIB per capita 2021': 'pib_2021'})

dummy_data = pd.read_csv('data/ips_brasil_municipios.csv')
dummy_data['Código IBGE'] = dummy_data['Código IBGE'].astype(str)

dummy_data['Município'] = dummy_data['Município'].apply(lambda x: str(unidecode(x)).upper()[:-5])
dummy_data['city_residency'] = dummy_data['UF'] + '_' + dummy_data['Município']

urban_pop = gpd.read_file('utils/urban/sidra_9923_PercPopResidUrbana_mun_22_polPolygon.shp')
urban_pop = urban_pop.set_index('Geocodigo').join(dummy_data.set_index('Código IBGE')[['UF']]).reset_index()

urban_pop['Nome'] = urban_pop['Nome'].apply(lambda x: str(unidecode(x)).upper())
urban_pop['city_residency'] = urban_pop['UF'] + '_' + urban_pop['Nome']

columns = ['PercPopRes', 'PopResidUr']
urban_pop = urban_pop[['city_residency'] + columns].rename(columns={'PercPopRes': 'percent_urban', 'PopResidUr': 'urban_pop'})

urban_area = pd.read_csv('data/urban_data_brazil.csv')

# join all different datasets
df = pd.concat([
    temp.set_index(['city_residency', 'date_first_symptoms']),
    prec.set_index(['city_residency', 'date_first_symptoms']),
    cases.set_index(['city_residency', 'date_first_symptoms'])[['dengue_inc', 'population', 'n_cases', 'country']],
], axis=1).reset_index()

df['state_residency'] = df['city_residency'].apply(lambda x: x.split('_')[0])
df['region'] = df['state_residency'].map(region_map)
df['year'] = df['date_first_symptoms'].dt.year
df['year'][df['country'] == 'Brazil'] = df[df['country'].isin(['Brazil', 'Peru'])]['date_first_symptoms'].apply(lambda x: x.year + 1 if x.month >= 8 else x.year)

df['year_childs'] = df['date_first_symptoms'].dt.year
df['year_childs'][df['country'] == 'Brazil'] = df[df['country'].isin(['Brazil', 'Peru'])]['date_first_symptoms'].apply(lambda x: x.year + 1 if x.month >= 8 else x.year)
df['month_childs'] = df['date_first_symptoms'].dt.month

df['year_childs'] = 'year_' + df['year_childs'].astype(str)
df['month_childs'] = 'month_' + df['month_childs'].astype(str)

df['month'] = df['date_first_symptoms'].dt.month
df = df[(df['year'] > 2001) & (df['year'] < 2025)]
df['dengue_inc'] = df['dengue_inc'].fillna(0)
df[(df['country'] == 'Peru') & (df['date_first_symptoms'] >= '2024-08-01')] = np.nan

df = df.dropna()
df = df.set_index(['year', 'city_residency']).join(immunity.set_index(['year', 'city_residency'])).reset_index()

df = df.set_index('city_residency').join(ips_data.set_index('city_residency')).join(urban_pop.set_index('city_residency')).reset_index()
df = df.set_index(['city_residency', 'year']).join(urban_area.set_index(['city_residency', 'year'])).reset_index()

df['month'] = 'month_' + df['month'].astype(int).astype(str) + df['region'].astype(str)
df['year'] = 'year_' + df['year'].astype(int).astype(str) + df['region'].astype(str)


df.to_csv('data/model_input_brazil_immunity_city.csv', index=False)