import pandas as pd
import numpy as np
import os


def calculate_incidence(temp_vars, temp_coefficients, data, coef_row):
    """Calculate incidence from coefficients for a specific bootstrap iteration"""
    # Get temperature data as matrix
    temp_matrix = data[temp_vars].values
    
    # Get coefficients as vector
    coef_vector = np.array([x for x in temp_coefficients.iloc[coef_row][temp_vars].values])
    
    # Calculate predictions using matrix multiplication
    log_incidence = temp_matrix @ coef_vector  # @ is matrix multiplication in Python
    pred_incidence = np.exp(log_incidence)
    
    # Add to data
    result = data.copy()
    result['pred_incidence'] = pred_incidence.astype(float)
    result['pred_cases'] = result['pred_incidence'] * result['population'] / 100000
    result['bootstrap_iteration'] = coef_row
    
    return result

def process_ensemble_no_save(temp_coefficients, temp_vars, ensemble_num, model):
    """Process a single ensemble without saving individual files"""
    print(f"Processing ensemble {ensemble_num}...")
    
    # Load ensemble-specific data
    ensemble_file = f"/gws/nopw/j04/cpdn_nonnerc/aaim/dengue/predict-all/ext-nat-ens{ensemble_num:03d}.csv"  # Change from NAT to ALL if necessary
    
    if not os.path.exists(ensemble_file):
        print(f"Warning: {ensemble_file} not found, skipping...")
        return None
    
    # Load ensemble data
    dengue_temp_ens = pd.read_csv(ensemble_file, parse_dates=['date_first_symptoms'])
    
    # Run all coefficient iterations
    n_boot = len(temp_coefficients)
    boot_results = []
    
    for i in range(n_boot):
        result = calculate_incidence(temp_vars, temp_coefficients, dengue_temp_ens, i)
        boot_results.append(result)
    
    # Combine all bootstrap results
    boot_results = pd.concat(boot_results, ignore_index=True)

    # Fitted
    dengue_temp = pd.read_csv("model_input_brazil_immunity_city.csv", parse_dates=['date_first_symptoms'])
    dengue_temp = dengue_temp[dengue_temp['date_first_symptoms'] >= '2023-08-01']
    boot_fitted = []

    for i in range(n_boot):
        result = calculate_incidence(temp_vars, temp_coefficients, dengue_temp, i)
        boot_fitted.append(result)
    
    # Combine all bootstrap results
    boot_fitted = pd.concat(boot_fitted, ignore_index=True)
    boot_fitted = boot_fitted.rename(
        columns={'pred_cases': 'fitted_cases', 'population': 'fitted_population'}
    )

    boot_results = boot_results.set_index([
        'region', 'city_residency', 'date_first_symptoms', 'bootstrap_iteration'
    ]).join(
        boot_fitted.set_index([
            'region', 'city_residency', 'date_first_symptoms', 'bootstrap_iteration'
        ])[['fitted_cases', 'fitted_population']]
    ).reset_index()

    # Aggregate by region and date
    final_results = (boot_results
                    .groupby(['region', 'date_first_symptoms', 'bootstrap_iteration'])
                    .apply(lambda x: pd.Series({
                        'predicted_incidence': x['pred_cases'].sum(),
                        'population': x['population'].sum(),
                        'fitted_incidence': x['fitted_cases'].sum(),
                        'fitted_population': x['fitted_population'].sum(),
                    }), include_groups=False)
                    .reset_index())

    final_results['model'] = model
    
    print(f"Ensemble {ensemble_num}: Completed {n_boot} coefficient iterations")
    
    return final_results


models = [
    {
        'model': 'lag1', 
        'filename': 'data/lag1_coef_state_blockboot1000.csv', 
        'vars': [
            "mean_2m_air_temp_degree1_lag1", "mean_2m_air_temp_degree2_lag1", "mean_2m_air_temp_degree3_lag1",
            "total_precipitation_lag1",
        ]
    },
    {
        'model': 'childs', 
        'filename': 'data/childs_coef_state_blockboot1000.csv', 
        'vars': [
            "mean_2m_air_temp_degree1_lag1", "mean_2m_air_temp_degree2_lag1", "mean_2m_air_temp_degree3_lag1",
            "mean_2m_air_temp_degree1_lag2", "mean_2m_air_temp_degree2_lag2", "mean_2m_air_temp_degree3_lag2",
            "mean_2m_air_temp_degree1_lag3", "mean_2m_air_temp_degree2_lag3", "mean_2m_air_temp_degree3_lag3",
            "total_precipitation_lag1", "total_precipitation_lag2", "total_precipitation_lag3", 
        ]
    },
    {
        'model': 'standard', 
        'filename': 'data/standard_coef_state_blockboot1000.csv', 
        'vars': [
            "mean_2m_air_temp_degree1_lag1", "mean_2m_air_temp_degree2_lag1", "mean_2m_air_temp_degree3_lag1",
            "mean_2m_air_temp_degree1_lag2", "mean_2m_air_temp_degree2_lag2", "mean_2m_air_temp_degree3_lag2",
            "mean_2m_air_temp_degree1_lag3", "mean_2m_air_temp_degree2_lag3", "mean_2m_air_temp_degree3_lag3",
            "total_precipitation_lag1", "total_precipitation_lag2", "total_precipitation_lag3", 
        ]
    },
    {
        'model': 'immunity', 
        'filename': 'data/with_immunity_coef_state_blockboot1000.csv', 
        'vars': [
            "mean_2m_air_temp_degree1_lag1", "mean_2m_air_temp_degree2_lag1", "mean_2m_air_temp_degree3_lag1",
            "mean_2m_air_temp_degree1_lag2", "mean_2m_air_temp_degree2_lag2", "mean_2m_air_temp_degree3_lag2",
            "mean_2m_air_temp_degree1_lag3", "mean_2m_air_temp_degree2_lag3", "mean_2m_air_temp_degree3_lag3",
            "total_precipitation_lag1", "total_precipitation_lag2", "total_precipitation_lag3", 
        ]
    },    
    {
        'model': 'lag45', 
        'filename': 'data/lag_45_coef_state_blockboot1000.csv', 
        'vars': [
            "mean_2m_air_temp_degree1_lag1", "mean_2m_air_temp_degree2_lag1", "mean_2m_air_temp_degree3_lag1",
            "mean_2m_air_temp_degree1_lag2", "mean_2m_air_temp_degree2_lag2", "mean_2m_air_temp_degree3_lag2",
            "mean_2m_air_temp_degree1_lag3", "mean_2m_air_temp_degree2_lag3", "mean_2m_air_temp_degree3_lag3",
            "mean_2m_air_temp_degree1_lag4", "mean_2m_air_temp_degree2_lag4", "mean_2m_air_temp_degree3_lag4",
            "mean_2m_air_temp_degree1_lag5", "mean_2m_air_temp_degree2_lag5", "mean_2m_air_temp_degree3_lag5",    
            "total_precipitation_lag1", "total_precipitation_lag2", "total_precipitation_lag3", 
            "total_precipitation_lag4", "total_precipitation_lag5"
        ]
    }, 
]

# Process all 525 ensembles and save to single file
print("Starting processing of 525 ensembles...")

all_results = []
results_summary = []
successful_ensembles = []
failed_ensembles = []

for model_iteration in models:
    print("Processing model: ", model_iteration['model'])

    # Load data
    temp_coefficients = pd.read_csv(model_iteration['filename']).iloc[:100] # Using just the first 100 rows for bootstrapping
    
    # Temperature variable names
    temp_vars = model_iteration['vars']

    for ensemble_num in range(525):  # 1 to 525
        try:
            result = process_ensemble_no_save(temp_coefficients, temp_vars, ensemble_num, model_iteration['model'])
            if result is not None:
                # Add ensemble number to results
                result['ensemble'] = ensemble_num
                all_results.append(result)
                
                successful_ensembles.append(ensemble_num)
                results_summary.append({
                    'ensemble': ensemble_num,
                    'status': 'success',
                    'n_predictions': len(result)
                })
            else:
                failed_ensembles.append(ensemble_num)
                results_summary.append({
                    'ensemble': ensemble_num,
                    'status': 'file_not_found',
                    'n_predictions': 0
                })
        except Exception as e:
            failed_ensembles.append(ensemble_num)
            results_summary.append({
                'ensemble': ensemble_num,
                'status': f'error: {str(e)}',
                'n_predictions': 0
            })
            print(f"Error processing ensemble {ensemble_num}: {e}")

# Combine all results and save to single file
if all_results:
    print("\nCombining all results into single master file...")
    master_results = pd.concat(all_results, ignore_index=True)
    
    # Save the combined results
    output_filename = f"/gws/nopw/j04/cpdn_nonnerc/aaim/dengue/sens_incidence_results_all_{len(successful_ensembles)}_nat_ensembles.csv"
    master_results.to_csv(output_filename, index=False)
    
    print(f"\n=== PROCESSING COMPLETE ===")
    print(f"Successfully processed: {len(successful_ensembles)} ensembles")
    print(f"Failed: {len(failed_ensembles)} ensembles") 
    print(f"Success rate: {len(successful_ensembles)/525*100:.1f}%")
    print(f"Total predictions: {len(master_results):,}")
    print(f"All results saved to: {output_filename}")
        
    if failed_ensembles:
        print(f"Failed ensembles: {failed_ensembles[:10]}{'...' if len(failed_ensembles) > 10 else ''}")
        
else:
    print("No ensembles were successfully processed!")
