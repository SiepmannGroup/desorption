    real::inv_dif_count

    inv_dif_count=real(new_count-last_count,dp)
    if (inv_dif_count.gt.0.5_dp) then
       inv_dif_count=1.0_dp/inv_dif_count
       ! The sequence of calculations is for avoiding overflow problems
       block_average=inv_dif_count*new_count*new_value-inv_dif_count*last_count*last_value
       last_value=new_value
       last_count=new_count
    end if
