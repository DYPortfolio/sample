function mixed_boundary_billiards_2022(varargin)
%Copyright Dmitry Yampolsky 2022

text_data_out_switch=false;
saveworkspace_switch=true;
savedata_switch = false;

mbb_main();

    function mbb_main()

        angle_alph = pi/6;
        a = 1;
        b = 0;
        maxrootN = 40;
        rootNs =  1:maxrootN; nKs = length(rootNs);
 
        ks = billiard_zeros(nKs,a,b);
      
        overlap_epsilon = 1e-9;
        angle_alph_side = tan(angle_alph);
        overlapF_inds=[];
        L=1;
        indexN=0;
        A = @(k)(1./ (((L./2)-sin(2.*k.*L)./(4.*k))) ^.5);

        %1st H
        hN = 1;
        kMax = ks(end);
        ijs = nchoosek(rootNs,2);
        nijs = length(ijs);

        %loop 1
        for ijctr = 1:nijs

            kctr_i = ijs(ijctr,1);
            kctr_j = ijs(ijctr,2);

            ki = (kctr_i);
            kj = (kctr_j);%overlapF1 input

            indexN=ijctr;
            kijref_ind(indexN,:)=[kctr_i,kctr_j];%indexes
            kijref(indexN,1:2) = [ks(kctr_i),ks(kctr_j)];%values

            ki_val = ks(kctr_i);
            kj_val = ks(kctr_j);

            ijHij_int(kctr_i,kctr_j) = ...
                (A(ki_val).^2 .* A(kj_val).^2 .*...
                ...
                ((ki_val.^2 + (kj_val).^2 .* cot(angle_alph).^2) .*...
                ...
                overlapF1(ki,ki,kj,kj,ki_val,ki_val,kj_val,kj_val) -...
                ...
                (ki_val.^2 + kj_val.^2).*csc(angle_alph).^2.*...
                ...
                overlapF1(ki,kj,ki,kj,ki_val,kj_val,ki_val,kj_val) + ...
                ...
                (kj_val.^2 + ki_val.^2.*cot(angle_alph).^2).*...
                ...
                overlapF1(kj,kj,ki,ki,kj_val,kj_val,ki_val,ki_val)))./2;

            indexed_ijHij_int2(indexN) =  ijHij_int(kctr_i,kctr_j) ;
        end% ij_ctr


        H_selection = (kijref(:,1).^2 + ((1./angle_alph_side).^2) .*(kijref(:,2).^2))...
            <  (max((1./angle_alph_side.^2).*(ks(end).^2),   ks(end).^2));

        kijref = kijref(H_selection,:);
        kijref_ind = kijref_ind(H_selection,:);

        indexed_ijHij_int2 = indexed_ijHij_int2(H_selection);
        [ijHij_int_sorted_cut, ijHij_int_sortOrder_cut] = sort(indexed_ijHij_int2);%sort ijHij

        %take root pairs at indexes ijHij_int_sortOrder_cut and their indexes in root list
        kijref_ind_sorted = kijref_ind(ijHij_int_sortOrder_cut,1:2);
        kijref_ind_sorted_val = kijref(ijHij_int_sortOrder_cut,1:2);

        hnnP = zeros(length(ijHij_int_sortOrder_cut),length(ijHij_int_sortOrder_cut));%init

        %2nd H
        hN = 2;
        size_kijref_ind_sorted=size(kijref_ind_sorted,1);
        size_kijref_ind_sorted_str = num2str(size(kijref_ind_sorted,1));%used to display progress

        overlapF1_handle = @overlapF1;

        size1 = size((kijref_ind_sorted),1);
        inds12(:,:,1)=repmat(1:size1,[size1,1]);
        inds12(:,:,2)=repmat(1:size1,[size1,1])';
        inds12 = reshape(inds12,size1^2,2);

        %loop 2
        parfor nctr0 = 1:size1^2

            nctr2 = inds12(nctr0,1);
            nctr1 = inds12(nctr0,2);

            ki = kijref_ind_sorted(nctr1,1);
            kj = kijref_ind_sorted(nctr1,2);
            ki_val = kijref_ind_sorted_val(nctr1,1);
            kj_val = kijref_ind_sorted_val(nctr1,2);
            kiP = kijref_ind_sorted(nctr2,1);
            kjP = kijref_ind_sorted(nctr2,2);
            kiP_val = kijref_ind_sorted_val(nctr2,1);
            kjP_val = kijref_ind_sorted_val(nctr2,2);

            hnnP(nctr0) =...
                (A(ki_val).*A(kiP_val).*A(kj_val).*A(kjP_val).*...
                ...
                ((kiP_val.^2 + kjP_val.^2.*cot(angle_alph).^2).*overlapF1_handle(ki,kiP,kj,kjP,ki_val,kiP_val,kj_val,kjP_val) - ...
                ...
                (kjP_val.^2 + kiP_val.^2.*cot(angle_alph).^2) .* overlapF1_handle(ki,kjP,kiP,kj,ki_val,kjP_val,kiP_val,kj_val) - ...
                ...
                (kiP_val.^2 + kjP_val.^2.*cot(angle_alph).^2) .* overlapF1_handle(kiP,kj,ki,kjP,kiP_val,kj_val,ki_val,kjP_val) + ...
                ...
                (kjP_val.^2 + kiP_val.^2.*cot(angle_alph).^2).*  overlapF1_handle(kj,kjP,ki,kiP,kj_val,kjP_val,ki_val,kiP_val)))./2  ;

        end%nctr1

        [h2_eig, h2_ev] = eig(hnnP);
        fprintf(" eigenvalues done\n")
        [sorted_eig, sorted_eigInds] = sort(diag(h2_ev));
        h2_eig = h2_eig(:,sorted_eigInds);

        output();

        
        function output()
            %...
        end
        
        clear%for multiple a/b's




        function overlapF1_out = overlapF1(kappaI1_ind, kappaII1_ind, kappaI2_ind, kappaII2_ind,kappaI1_val, kappaII1_val, kappaI2_val, kappaII2_val)
            %function input are indexes not values, overlapF2 refers to kijref_ind_sorted, returns values

            overlapF_inds = [kappaI1_ind kappaII1_ind kappaI2_ind kappaII2_ind];%for each overlapF2 input below, indexes are same, only sign changes,
            %which gets delivered through input rather than global

            kappaI1 = kappaI1_val;
            kappaII1 = kappaII1_val;
            kappaI2 = kappaI2_val;
            kappaII2 = kappaII2_val;

            overlapF1_out = (1./8) .* (...
                + overlapF2(+kappaI1, +kappaII1, +kappaI2, +kappaII2)...
                - overlapF2(+kappaI1, +kappaII1, -kappaI2, +kappaII2)   ...
                + overlapF2(+kappaI1, -kappaII1, +kappaI2, -kappaII2)...
                - overlapF2(+kappaI1, -kappaII1, -kappaI2, -kappaII2)   ...
                - overlapF2(+kappaI1, +kappaII1, +kappaI2, -kappaII2)...
                + overlapF2(+kappaI1, +kappaII1, -kappaI2, -kappaII2)   ...
                - overlapF2(+kappaI1, -kappaII1, +kappaI2, +kappaII2)...
                + overlapF2(+kappaI1, -kappaII1, -kappaI2, +kappaII2)) ;

        end

        function  overlapF2_out = overlapF2(qI1, qII1, qI2, qII2)

            overlapF2_condition = 0;
            overlapF2_out=0;

            if (abs(qI1 + qI2 + qII1 + qII2) > overlap_epsilon * kMax) && (abs(qI1 + qII1) > overlap_epsilon * kMax ) && (abs(qI2 + qII2) > overlap_epsilon * kMax)

                overlapF2_condition = 1;
                overlapF2_out = ...
                    -(qI2./((qI1 + qII1).*(qI2 + qII2).*(qI1 + qI2 + qII1 + qII2))) - qII2./((qI1 + qII1).*(qI2 + qII2).*(qI1 + qI2 + qII1 + qII2)) + ...
                    cos(qI1 + qII1)./((qI1 + qII1).*(qI2 + qII2)) - cos(qI1 + qI2 + qII1 + qII2)./((qI2 + qII2).*(qI1 + qI2 + qII1 + qII2));

            elseif (abs(qI1 + qII1 )<= overlap_epsilon * kMax ) && (abs(qI2 + qII2) > overlap_epsilon * kMax )

                overlapF2_condition = 2;
                overlapF2_out= (1 - cos(qI2 + qII2))./(qI2 + qII2).^2;

            elseif (abs(qI1 + qII1 ) > overlap_epsilon * kMax ) && (abs(qI2 + qII2) <= overlap_epsilon * kMax )

                overlapF2_condition = 3;
                overlapF2_out = (-1 + cos(qI1 + qII1) + (qI1 + qII1).*sin(qI1 + qII1))./(qI1 + qII1).^2;

            elseif (abs(qI1 + qII1 ) <= overlap_epsilon * kMax ) && (abs(qI2 + qII2) <= overlap_epsilon * kMax )

                overlapF2_condition = 4;
                overlapF2_out = 1/2;

            elseif  (abs(qI1 + qI2 + qII1 + qII2)  <= (overlap_epsilon * kMax)) && (abs(qI2 + qII2) > overlap_epsilon * kMax)

                overlapF2_condition = 5;
                overlapF2_out = (1 - cos(qI2 + qII2))./(qI2 + qII2).^2;
            else
                overlapF2_condition = 6;
                error('valid condition in overlapF2_out not found');

            end
        end
    end

end
