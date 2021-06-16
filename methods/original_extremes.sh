pathwayIN='../../../../australia_climate/prec'
CMIP6='SSP245_r1i1p1f1_K_1850_2100.nc'

for model in CanESM5; do
    cdo -L -mulc,86400 -selyear,1851/2100 -mul mask_${model}.nc \
        ${pathwayIN}/prec_${model}_${CMIP6} \
        ${model}/prec_${model}_1851-2100.nc
        
    ### Calculate percentiles
    cdo timmin ${model}/prec_${model}_1851-2100.nc \
               ${model}/prec_minfile.nc
    cdo timmax ${model}/prec_${model}_1851-2100.nc \
               ${model}/prec_maxfile.nc
    cdo timpctl,5 ${model}/prec_${model}_1851-2100.nc \
        ${model}/prec_minfile.nc \
        ${model}/prec_maxfile.nc \
        ${model}/prec_5_percentile.nc
    cdo timpctl,95 ${model}/prec_${model}_1851-2100.nc \
        ${model}/prec_minfile.nc \
        ${model}/prec_maxfile.nc \
        ${model}/prec_95_percentile.nc

    ### Make percentile masks
    cdo -L -lec,0 -sub \
        ${model}/prec_${model}_1851-2100.nc \
        ${model}/prec_5_percentile.nc \
        ${model}/prec_mask_only_5.nc
    cdo -L -gec,0 -sub \
        ${model}/prec_${model}_1851-2100.nc \
        ${model}/prec_95_percentile.nc \
        ${model}/prec_mask_only_95.nc
    cdo -L -mulc,-1 -subc,1 ${model}/prec_mask_only_5.nc \
        ${model}/prec_mask_except_5.nc
    cdo -L -mulc,-1 -subc,1 ${model}/prec_mask_only_95.nc \
        ${model}/prec_mask_except_95.nc

    ### Get 5/95 percentile in original data        
    cdo mul ${model}/prec_${model}_1851-2100.nc \
        ${model}/prec_mask_only_5.nc \
        ${model}/prec_original_5.nc
    cdo mul ${model}/prec_${model}_1851-2100.nc \
        ${model}/prec_mask_only_95.nc \
        ${model}/prec_original_95.nc  
        
    for method in CDFt; do    
        ### Set 5/95 percentile to zero in corrected data      
        cdo mul ${method}/${model}/${method}_prec_${model}_COR_1851-2100.nc \
            ${model}/prec_mask_except_5.nc \
            ${method}/${model}/prec_bc_except_5.nc
        cdo mul ${method}/${model}/${method}_prec_${model}_COR_1851-2100.nc \
            ${method}/${model}/prec_mask_except_95.nc \
            ${method}/${model}/prec_bc_except_95.nc

        ### Corrected data with original 5 percentile data
        cdo add ${model}/prec_original_5.nc \
            ${method}/${model}/prec_bc_except_5.nc \
            ${method}/${model}/prec_no_low_extremes.nc
            
        ### Corrected data with original 95 percentile data
        cdo add ${model}/prec_original_95.nc \
            ${method}/${model}/prec_bc_except_95.nc \
            ${method}/${model}/prec_no_high_extremes.nc

        ### Corrected data with original 5 and 95 percentile data    
        cdo mul ${method}/${model}/prec_no_low_extremes.nc \
            ${model}/prec_mask_except_95.nc \
            ${method}/${model}/prec_no_low_except_95.nc
        cdo add ${model}/prec_original_95.nc \
            ${method}/${model}/prec_no_low_except_95.nc \
            ${method}/${model}/prec_no_extremes.nc
            
        done
done
