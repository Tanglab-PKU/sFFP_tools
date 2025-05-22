import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.preprocessing import StandardScaler
import numpy as np
import joblib  # For model saving
from sklearn import tree  # For decision tree visualization

#%% Data Loading
# Define file paths for experimental datasets
data_paths = [
    "D:/MSdata/231215_10MIX/250304Search/ADK-SDA(DESTHY)/reports/ADK-SDA(DESTHY)_con_2025.03.05_XLs_Jnoted_Dis_FFPiso-discrete-noted_FDR.xlsx",
    "D:/MSdata/231215_10MIX/250304Search/NSP5-SDA(DESTHY)/reports/NSP5-SDA(DESTHY)_con_2025.03.05_XLs_Jnoted_Dis_FFPiso-discrete-noted_FDR.xlsx",
    "D:/MSdata/231215_10MIX/250304Search/CA-SDA(DESTHY)/reports/CA-SDA(DESTHY)_con_2025.03.05_XLs_Jnoted_Dis_FFPiso-discrete-noted_FDR.xlsx",
    "D:/MSdata/231215_10MIX/250304Search/LF-SDA(DESTHY)/reports/LF-SDA(DESTHY)_con_2025.03.05_XLs_Jnoted_Dis_FFPiso-discrete-noted_FDR.xlsx",
    "D:/MSdata/231215_10MIX/250304Search/GFP-SDA(DESTHY)/reports/GFP-SDA(DESTHY)_con_2025.03.05_XLs_Jnoted_Dis_FFPiso-discrete-noted_FDR.xlsx",
    "D:/MSdata/231215_10MIX/250304Search/BSA-SDA(DESTHY)/reports/BSA-SDA(DESTHY)_con_2025.03.05_XLs_Jnoted_Dis_FFPiso-discrete-noted_FDR.xlsx",
    # ... (other paths remain same)
]

# Load and concatenate data with stratified sampling
dataset = pd.concat(
    [pd.read_excel(path).iloc[:2500] for path in data_paths], 
    axis=0
).fillna(0)
print(f"Initial dataset shape: {dataset.shape}")

#%% Feature Engineering
# Calculate peptide spectral matches (PSMs)
dataset["Peptide_PSMs"] = dataset.groupby("Peptide")["Peptide"].transform("size")

# Feature selection (active features marked with *)
feature_columns = [
    'SD_A_base_alpha','SD_B_base_alpha',    
    'Ion_Count_A_base_alpha','Ion_Count_B_base_alpha', #  0.7900
    
    # 'sum_A_y0_α_RelInt','sum_B_y0_α_RelInt', # 0.7350(FDR-good!(not better than independently use))
    # 'A_y0_pepB_RelInt_base_alpha','B_y0_pepA_RelInt_base_alpha', # 0.5100
    'A_y0_α_RelInt_base_alpha',
    'B_y0_ω_RelInt_base_alpha', # 0.7175(FDR-good!)
    'A_y0_ω_RelInt_base_alpha',
    'B_y0_α_RelInt_base_alpha', #  0.6625(FDR-good!)
    # 'Total_Relative_Intensity_A_alpha',
    # 'Total_Relative_Intensity_B_sz',
    # 'Total_Relative_Intensity_A_sz',
    # 'Total_Relative_Intensity_B_alpha',
    # "Total_Relative_Intensity_B_counter_N_all","Total_Relative_Intensity_A_counter_N_all"

    "Peptide_PSMs",
    # 'Average_Relative_Intensity_base','Average_Relative_Intensity_alpha_sz',
    # 'Q-value','Q-value_CSM',
    # 'Score', # 0.8300(FDR-BAD!! BUT GOOD Accuracy)
    # 'Charge',
    # 'isFilterIn','isComplexSatisfied',       'isRetainHighLevel',
    # 'Target_Decoy', #Warn!
    # ... (other features)
]

X = dataset[feature_columns].fillna(0)
y = dataset['waste'].fillna(0)

#%% Data Partitioning
train_tags = [            
            'ADK_L',
            'GFP_L',
            'CA_L',
            'LF_L',
            'NSP5_L',]  # Training set identifiers
test_tags = ['BSA_L']            # Test set identifier

# Create boolean masks using regular expressions
train_mask = dataset['Title'].str.contains('|'.join(train_tags), case=False, na=False)
test_mask = dataset['Title'].str.contains('|'.join(test_tags), case=False, na=False)

# Ensure mutual exclusivity
train_mask = ~test_mask & train_mask

# Split datasets
X_train, X_test = X[train_mask], X[test_mask]
y_train, y_test = y[train_mask], y[test_mask]

#%% Model Configuration
param_grid = {
    'n_estimators': [500],  # Optimal number of trees
    'max_depth': [4],       # Controls model complexity
}

# Initialize grid search with 3-fold CV
grid_search = GridSearchCV(
    estimator=RandomForestClassifier(random_state=42),
    param_grid=param_grid,
    cv=3,
    n_jobs=1,
    verbose=2
)

#%% Model Training
grid_search.fit(X_train, y_train)
best_model = grid_search.best_estimator_

# Save trained model
joblib.dump(best_model, r'D:\MSdata\PRIDE_UPLOAD-202505\Related Code\SDA-random_forest_model.pkl')  # you can load model by loaded_model = joblib.load('SDA-random_forest_model.pkl')

#%% Model Evaluation
y_pred = best_model.predict(X_test)
print(f"Best Parameters: {grid_search.best_params_}")
print(f"Test Accuracy: {accuracy_score(y_test, y_pred):.4f}")

# Feature importance analysis
importances = pd.DataFrame({
    'Feature': X.columns,
    'Importance': best_model.feature_importances_
}).sort_values('Importance', ascending=False)


#%% Results Export
predictions = dataset[test_mask].copy()
predictions['Predicted_Class'] = y_pred
predictions['Prediction_Probability'] = best_model.predict_proba(X_test)[:,1]

output_path = r"D:\MSdata\PRIDE_UPLOAD-202505\Related Code\BSA_L_Predictions.xlsx"
predictions.to_excel(output_path, index=False)
