#Model Optimizer -------------------------------------------------------------------------------------------------------

def create_OP_Model(trial):

    #To save memory
    clear_session()
    
    # Hidden Layers Optimization Range Settings: 

    #Image layers suggestions
    n_image_layers = trial.suggest_int('N_Image_layers', 1, 10)

    #Energy layers suggestions
    #n_en_layers = trial.suggest_int('N_Energy_layers', 1, 3)

    #Radius Layers suggestions
    #n_rad_layers = trial.suggest_int('N_Radius_layers', 1, 3)

    #Before Dropout Layers (BDO)
    n_bdo_layers = trial.suggest_int('N_bdo_layers', 1, 10)

    #After Dropout (ADO) Layers
    n_ado_layers = trial.suggest_int('N_ado_layers', 1, 10)

    #Adam Optimizer Hyperparams: ----------------------------------------------------------------
    #Learning Rate
    lr_rate = trial.suggest_float('learning_rate', 1e-4, 1e-3)
    #Momentum estimators beta_1 (mean) and beta_2 (variance)
    #beta_1 = trial.suggest_float('momentum_beta_1', 0.8, 0.95)
    #beta_2 = trial.suggest_float('momentum_beta_2', 0.98, 0.999)

    #Dropout rate: ------------------------------------------------------------------------------
    do_rate = trial.suggest_float('Drop_out_rate', 0, 0.7)

    #lists with the number of neurons to be saved in each Optimized Layer
    n_units_image_ly = []
    #n_units_en_ly = []
    #n_units_rad_ly = []
    n_units_bdo_ly = []
    n_units_ado_ly = []

#-----------------------------------Architecture and Nº Neurons Suggestion------------------------
    
    #Input Layers --------------------------------------------------------------------------------
    # Image Input Layer
    image_dim = (10,10)
    input_im_shape = image_dim[0]*image_dim[1]
    input_image = Input(shape = (input_im_shape,), name = 'image_input')

    #Energy Input Layer
    input_energy = Input(shape = (1,), name = 'energy_input')

    #Position Input Layer
    input_position = Input(shape = (1,), name = 'position_input')

    #Hidden Layers for pixels values with decreasing number of neurons: -------------------------
    n_neurons_image_first_ly = trial.suggest_int('n_neurons_1st_image_hidden_ly', 80, 95) #Nº of maximum neurons in the first hidden layer of the input image pixel values
    n_neurons_image_last_ly = trial.suggest_int('n_neurons_last_image_hidden_ly', 20, 30)  #Nº of minimum neurons in the last hidden layer of the input image pixel values
    for i in range(n_image_layers):
        if i == 0:
            n_imag_neurons = trial.suggest_int(f'n_neurons_im_ly_{i}', n_neurons_image_last_ly, n_neurons_image_first_ly) #first Image Hidden Layer with a maximum of
            image_ly = Dense(n_imag_neurons, activation = 'gelu', name = f'im_ly_{i}')(input_image)
        else:
            n_imag_neurons = trial.suggest_int(f'n_neurons_im_ly_{i}', n_neurons_image_last_ly, n_units_image_ly[i-1])  #Garantiza que el numero de neuronas disminuye, o es el mismo, con un minimo de n_neurons_image_last_ly!
            image_ly = Dense(n_imag_neurons, activation = 'gelu', name = f'im_ly_{i}')(image_ly)
        n_units_image_ly.append(n_imag_neurons)

    #Concatenate layer --------------------------------------------------------------------------
    x = concatenate([image_ly, input_energy, input_position])
    
    # Before-Drop-Out layers --------------------------------------------------------------------
    n_neurons_bdo_last_ly = trial.suggest_int('n_neurons_last_bdo_hidden_ly', n_units_image_ly[-1] - 10, n_units_image_ly[-1])  #Nº of minimum neurons in the last BDO hidden layer as function of the neurons in the image_last_ly
    for i in range(n_bdo_layers):
        if i == 0:
            n_neurons = trial.suggest_int(f'n_neurons_bdo_ly_{i}', n_neurons_bdo_last_ly, n_units_image_ly[-1] + 2)  #Nº max neurons = neurons_last_imag_ly + 2 (energy + position)
            x = Dense(n_neurons, activation = 'gelu', name = f'bdo_layer{i}')(x) # 1st hidden layer after concatenate: Maximum neurons =  min(imag_units) of last layer + 2 units (energy and position)
        else:
            n_neurons = trial.suggest_int(f'n_neurons_bdo_ly_{i}', n_neurons_bdo_last_ly, n_units_bdo_ly[i-1])
            x = Dense(n_neurons, activation = 'gelu', name = f'bdo_layer{i}')(x)
        n_units_bdo_ly.append(n_neurons)
        

    # Dropout Layer ---------------------------------------------------------
    x = Dropout(rate = do_rate)(x)

    # After-Drop-Out layers --------------------------------------------------
    n_neurons_ado_last_ly = trial.suggest_int('n_neurons_last_ado_hidden_ly', 4, n_units_bdo_ly[-1])  #Nº of minimum neurons in the last ADO hidden layer as function of the neurons in the bdo_last_ly
    for i in range(n_ado_layers):
        if i == 0:
            n_neurons = trial.suggest_int(f'n_neurons_ado_ly_{i}', n_neurons_ado_last_ly, n_neurons)
        else:
            n_neurons = trial.suggest_int(f'n_neurons_ado_ly_{i}', n_neurons_ado_last_ly, n_units_ado_ly[i-1])
        n_units_ado_ly.append(n_neurons)
        x = Dense(n_neurons, activation = 'gelu', name = f'ado_ly_{i}')(x)

    # Output Layer -----------------------------------------------------------------------    
    output_layer = Dense(2, activation = 'sigmoid', name = 'nu_predict')(x)
    
    #--------------------------------------------------------------------------------------

    model = Model(inputs = [input_image, input_energy, input_position], outputs = output_layer)
    
    model.compile(optimizer = tf.keras.optimizers.Adam(learning_rate = lr_rate),
                  loss = tf.keras.losses.BinaryCrossentropy(),
                  metrics = ['AUC', 'accuracy'])

    return model

    #--------------------------------------------------------------------------------------

def objective(trial):

    model = create_OP_Model(trial)
    
    history = model.fit(x = {'image_input': pixel_train1_transf, 'energy_input': energy_train1_transf, 'position_input': position_train1_transf}, 
                    y = {'nu_predict':labels_train1}, 
                    epochs = 200, 
                    batch_size = 2000,
                    validation_data=([pixel_val_transf, energy_val_transf, position_val_transf], labels_val),
                    callbacks=[tf.keras.callbacks.EarlyStopping(patience = 15, min_delta = 1e-3, monitor="val_loss")],
                    shuffle = True)
    
    X_sig_pred = model.predict([pixel_test1_transf_sig, energy_test1_transf_sig, position_test1_transf_sig])
    X_bkg_pred = model.predict([pixel_test1_transf_bkg, energy_test1_transf_bkg, position_test1_transf_bkg])

    # Criterio de Optimización: 1) Disminuir error de predicción (BCE), 2) Aumentar AUC (probar)
    #1)BCE
    #BCE = tf.keras.losses.BinaryCrossentropy()

    #2)AUC
    x_pred = np.concatenate((X_sig_pred, X_bkg_pred))             # Predicted
    x_exp = np.concatenate((labels_test1_sig, labels_test1_bkg))  # Expected

    # Classification overall Model
    fpr, tpr, thresholds = roc_curve(np.concatenate((x_exp[:,0], x_exp[:,1])), np.concatenate((x_pred[:,0], x_pred[:,1])))
    auc_val = auc(fpr, tpr)

    return auc_val

# Initialize Optimizer -----------------------------------------------------------------
study = optuna.create_study(directions = ["maximize"])
study.optimize(objective, n_trials = 50)
best_trial = study.best_trial

print("Best Hyperparms.:", study.best_trial.params) #To see the best 
print("Value: ", best_trial.value)

# Run Best Model ------------------------------------------------------------------------
best_model = create_OP_Model(best_trial)
history = best_model.fit(x = {'image_input': pixel_train1_transf, 'energy_input': energy_train1_transf, 'position_input': position_train1_transf},
                    y = {'nu_predict':labels_train1}, 
                    epochs = 200, 
                    batch_size = 2000,
                    validation_data=([pixel_val_transf, energy_val_transf, position_val_transf], labels_val),
                    callbacks=[tf.keras.callbacks.EarlyStopping(patience = 15, min_delta = 1e-3, monitor="val_loss")],
                    shuffle = True)
    
#Plots of Train and Validation -----------------------------------------------------------
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
fig, axes = plt.subplots(1,2, figsize=(15, 5))

axes[0].plot(history.history['loss'], 'b:.', label = ' Train Loss')
axes[0].plot(history.history['val_loss'], 'r:.', label = ' Val. Loss')
axes[0].plot(history.history['accuracy'], 'g:.',  label = ' Train Accuracy')
axes[0].plot(history.history['val_accuracy'], 'c:.', label = ' Val. Accuracy')
axes[1].plot(history.history['AUC'], 'b:.',  label = ' Train AUC')
axes[1].plot(history.history['val_AUC'], 'r:.', label = ' Val. AUC')

axes[0].set_title('Loss and Accuracy Curves', fontsize = 10)
axes[0].set_xlabel('epoch')
axes[0].legend(loc = 'best')
axes[1].set_title('AUC Value Curve', fontsize = 10)
axes[1].set_xlabel('epoch')
axes[1].legend(loc = 'lower right')

#axes[0].set_ylim(0, 1.2)
#axes[1].set_ylim(0, 1.05)

#axes[0].set_xlim(-2, len(history.history['loss']))
#axes[1].set_xlim(-2, len(history.history['loss']))

fig.suptitle(f'Train and Validation Progress - Energy : [2.5, {energy_cut_train}] MeV - R < 5500 (mm)', backgroundcolor='blue', color='white', fontsize = 10)
#plt.savefig('figs/Optimized Models/10x10/Train_Val_progress.png', format = 'png', bbox_inches = 'tight')
plt.show()

#Observation of model Predictions ------------------------------------------------------------
X_sig_pred = best_model.predict([pixel_test2_transf_sig, energy_test2_transf_sig, position_test2_transf_sig])
X_bkg_pred = best_model.predict([pixel_test2_transf_bkg, energy_test2_transf_bkg, position_test2_transf_bkg])

sn.reset_orig
bins = 65
fig, axes = plt.subplots(1, 2,  figsize=(22, 7))

sn.histplot(X_sig_pred[:,0], bins = bins, label = 'SAS', color = 'r', alpha = 0.6, ax = axes[0], log = True)
sn.histplot(X_sig_pred[:,1], bins = bins, label = 'SAB', color = 'b', alpha = 0.6, ax = axes[0])
sn.histplot(X_bkg_pred[:,0], bins = bins, label = 'BAS',color = 'b', alpha = 0.6, ax = axes[1], log = True)
sn.histplot(X_bkg_pred[:,1], bins = bins, label = 'BAB', color = 'r', alpha = 0.6, ax = axes[1])
axes[0].set_title('Probability Components for Solar $ν$ Predictions', fontsize = 10)
axes[1].set_title('Probability Components for Background Predictions', fontsize = 10)
axes[0].legend(loc = 'best')
axes[1].legend(loc = 'best')
plt.suptitle(f'Model Predictions for Solar $^8$B-$ν_e$ events and Background events - Energy : [2.5, {energy_cut_train}] MeV - R < 5500 (mm)', backgroundcolor='blue', color='white')
#plt.savefig('figs/Optimized Models/10x10/probability_predictions1.png', format = 'png', bbox_inches = 'tight')
plt.show()

#ROC Curve ---------------------------------------------------------------------------------
X_pred = np.concatenate((X_sig_pred, X_bkg_pred))              #Predicted
X_exp = np.concatenate((labels_test2_sig, labels_test2_bkg))   #Expected

# Classification for solar_nu events
fpr_sig, tpr_sig, thresholds_sig = roc_curve(X_exp[:,0], X_pred[:,0])
auc_sig = auc(fpr_sig, tpr_sig)

# Classification for bkg events
fpr_bkg, tpr_bkg, thresholds_bkg = roc_curve(X_exp[:,1], X_pred[:,1])
auc_bkg = auc(fpr_bkg, tpr_bkg)

# Classification overall Model
fpr, tpr, thresholds = roc_curve(np.concatenate((X_exp[:,0], X_exp[:,1])), np.concatenate((X_pred[:,0], X_pred[:,1])))
auc_val = auc(fpr, tpr)
fig, axes = plt.subplots(1, 2,  figsize=(14, 5))


axes[0].plot([0, 1], [0, 1], 'k--') #x=y
axes[0].plot(fpr_sig, tpr_sig, 'r:,', label = 'sig ROC (area = {:.4f})'.format(auc_sig))
axes[0].plot(fpr_bkg, tpr_bkg, 'b:,', label = 'bkg ROC (area = {:.4f})'.format(auc_bkg))

axes[1].plot([0, 1], [0, 1], 'k--') #x=y
axes[1].plot(fpr, tpr, 'b:,', label = 'ROC over all Model (area = {:.4f})'.format(auc_val))

axes[0].set_ylabel('True Positive Rate')
axes[0].set_xlabel('False Positive Rate')
axes[1].set_ylabel('True Positive Rate')
axes[1].set_xlabel('False Positive Rate')

axes[0].legend(loc = 'lower right')
axes[1].legend(loc = 'lower right')

fig.suptitle(f'ROC curve for Model Predictions - Energy : [2.5, {energy_cut_train}] MeV - R < 5500 (mm)', fontsize = 10)

#plt.savefig('figs/Optimized Models/10x10/ROC.png', format = 'png', bbox_inches = 'tight')
plt.show()