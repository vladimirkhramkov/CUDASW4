#ifndef MANYPASS_FLOAT_KERNEL_CUH
#define MANYPASS_FLOAT_KERNEL_CUH

#include <cuda_fp16.h>

template <int numRegs, class PositionsIterator> 
struct ManyPassFloat{
    using BLOSUM62_SMEM = float[21][21];

    static constexpr int group_size = 32;

    static constexpr float negInftyFloat = -10000.0f;

    BLOSUM62_SMEM& shared_BLOSUM62;

    int numSelected;
    int length_2;
    float gap_open;
    float gap_extend;
    PositionsIterator d_positions_of_selected_lengths;
    const char* devChars;
    short2* devTempHcol2;
    short2* devTempEcol2;
    const size_t* devOffsets;
    const size_t* devLengths;

    __device__
    ManyPassFloat(
        BLOSUM62_SMEM& shared_BLOSUM62_,
        const char* devChars_,
        short2* devTempHcol2_,
        short2* devTempEcol2_,
        const size_t* devOffsets_,
        const size_t* devLengths_,
        PositionsIterator d_positions_of_selected_lengths_,
        int length_2_,
        float gap_open_,
        float gap_extend_
    ) : shared_BLOSUM62(shared_BLOSUM62_),
        devChars(devChars_),
        devTempHcol2(devTempHcol2_),
        devTempEcol2(devTempEcol2_),
        devOffsets(devOffsets_),
        devLengths(devLengths_),
        d_positions_of_selected_lengths(d_positions_of_selected_lengths_),
        length_2(length_2_),
        gap_open(gap_open_),
        gap_extend(gap_extend_)
    {



    }


    __device__
    void checkHEindex(int x, int line) const{
        // if(x < 0){printf("line %d\n", line);}
        // assert(x >= 0); //positive index
        // assert(2*(blockDim.x/group_size)*blockIdx.x * length_2 <= base_3 + x);
        // assert(base_3+x < 2*(blockDim.x/group_size)*(blockIdx.x+1) * length_2);
    };

    __device__
    void initial_calc32_local_affine_float(const int value, char query_letter, float& E, float& penalty_here31, float penalty_diag, float penalty_left, float& maximum, 
        int (&subject)[numRegs], 
        float (&penalty_here_array)[numRegs],
        float (&F_here_array)[numRegs]
    ) const{
        const float* const sbt_row = shared_BLOSUM62[query_letter];

        const float score2_0 = sbt_row[subject[0]];
        //score2.y = sbt_row[subject1[0].x];
        float penalty_temp0 = penalty_here_array[0];
        if (!value || (threadIdx.x%group_size)) E = max(E+gap_extend, penalty_left+gap_open);
        F_here_array[0] = max(F_here_array[0]+gap_extend, penalty_here_array[0]+gap_open);
        const float temp0_0 = max(max(penalty_diag + score2_0, max(E, F_here_array[0])), 0.f);
        penalty_here_array[0] = temp0_0;
        maximum = max(temp0_0, maximum);

        const float score2_1 = sbt_row[subject[1]];
        //score2.y = sbt_row[subject1[0].y];
        float penalty_temp1 = penalty_here_array[1];
        E = max(E+gap_extend, penalty_here_array[0]+gap_open);
        F_here_array[1] = max(F_here_array[1]+gap_extend, penalty_here_array[1]+gap_open);
        penalty_here_array[1] = max(max(penalty_temp0 + score2_1, max(E, F_here_array[1])), 0.f);
        const float temp0_1 = penalty_here_array[1];
        maximum = max(temp0_1,maximum);

        #pragma unroll
        for (int i=1; i<numRegs/2; i++) {
            const float score2_2i = sbt_row[subject[2*i]];
            //score2.y = sbt_row[subject1[i].x];
            penalty_temp0 = penalty_here_array[2*i];
            E = max(E+gap_extend, penalty_here_array[2*i-1]+gap_open);
            F_here_array[2*i] = max(F_here_array[2*i]+gap_extend, penalty_here_array[2*i]+gap_open);
            penalty_here_array[2*i] = penalty_here_array[2*i] = max(max(penalty_temp1 + score2_2i, max(E, F_here_array[2*i])), 0.f);
            const float temp0_2 = penalty_here_array[2*i]; 
            maximum = max(temp0_2,maximum);

            const float score2_2i1 = sbt_row[subject[2*i+1]];
            //score2.y = sbt_row[subject1[i].y];
            penalty_temp1 = penalty_here_array[2*i+1];
            E = max(E+gap_extend, penalty_here_array[2*i]+gap_open);
            F_here_array[2*i+1] = max(F_here_array[2*i+1]+gap_extend, penalty_here_array[2*i+1]+gap_open);
            penalty_here_array[2*i+1] = penalty_here_array[2*i+1] = max(max(penalty_temp0 + score2_2i1, max(E, F_here_array[2*i+1])), 0.f);
            const float temp0_3 = penalty_here_array[2*i+1];
            maximum = max(temp0_3,maximum);
        }

        penalty_here31 = penalty_here_array[numRegs-1];
        E = max(E+gap_extend, penalty_here31+gap_open);
        #pragma unroll //UNROLLHERE
        for (int i=0; i<numRegs; i++) F_here_array[i] = max(F_here_array[i]+gap_extend, penalty_here_array[i]+gap_open);
    };

    __device__
    void calc32_local_affine_float(char query_letter, float& E, float& penalty_here31, float penalty_diag, float& maximum, 
        int (&subject)[numRegs],
        float (&penalty_here_array)[numRegs],
        float (&F_here_array)[numRegs]
    ) const{
        const float* const sbt_row = shared_BLOSUM62[query_letter];

        const float score2_0 = sbt_row[subject[0]];
        float penalty_temp0 = penalty_here_array[0];
        penalty_here_array[0] = max(max(penalty_diag + score2_0, max(E, F_here_array[0])), 0.f);
        //if (penalty_here_array[0] > maximum) maximum = penalty_here_array[0];
        float penalty_temp1 = penalty_here_array[0]+gap_open;
        E = max(E+gap_extend, penalty_temp1);
        F_here_array[0] = max(F_here_array[0]+gap_extend, penalty_temp1);

        const float score2_1 = sbt_row[subject[1]];
        penalty_temp1 = penalty_here_array[1];
        penalty_here_array[1] = max(max(penalty_temp0 + score2_1, max(E, F_here_array[1])), 0.f);
        //if (penalty_here_array[1]  > maximum) maximum = penalty_here_array[1] ;
        penalty_temp0 = penalty_here_array[1]+gap_open;
        E = max(E+gap_extend, penalty_temp0);
        F_here_array[1] = max(F_here_array[1]+gap_extend, penalty_temp0);
		maximum = max(maximum,max(penalty_here_array[0],penalty_here_array[1]));


        #pragma unroll
        for (int i=1; i<numRegs/2; i++) {
            const float score2_2i = sbt_row[subject[2*i]];
            penalty_temp0 = penalty_here_array[2*i];
            penalty_here_array[2*i] = max(max(penalty_temp1 + score2_2i, max(E, F_here_array[2*i])), 0.f);
            //if (penalty_here_array[2*i] > maximum) maximum = penalty_here_array[2*i];
            penalty_temp1 = penalty_here_array[2*i]+gap_open;
            E = max(E+gap_extend, penalty_temp1);
            F_here_array[2*i] = max(F_here_array[2*i]+gap_extend, penalty_temp1);

            const float score2_2i1 = sbt_row[subject[2*i+1]];
            penalty_temp1 = penalty_here_array[2*i+1];
            penalty_here_array[2*i+1] = max(max(penalty_temp0 + score2_2i1, max(E, F_here_array[2*i+1])), 0.f);
            //if (penalty_here_array[2*i+1] > maximum) maximum = penalty_here_array[2*i+1];
            penalty_temp0 = penalty_here_array[2*i+1]+gap_open;
            E = max(E+gap_extend, penalty_temp0);
            F_here_array[2*i+1] = max(F_here_array[2*i+1]+gap_extend, penalty_temp0);
			maximum = max(maximum,max(penalty_here_array[2*i],penalty_here_array[2*i+1]));

        }

		//for (int i=0; i<numRegs/4; i++)
		//	maximum = max(maximum,max(max(penalty_here_array[4*i],penalty_here_array[4*i+1]), max(penalty_here_array[4*i+2],penalty_here_array[4*i+3])));
        penalty_here31 = penalty_here_array[numRegs-1];
    }

    __device__    
    void init_penalties_local(int value, float& penalty_diag, float& penalty_left, 
        float (&penalty_here_array)[numRegs], 
        float (&F_here_array)[numRegs]
    ) const{
        penalty_left = negInftyFloat;
        penalty_diag = negInftyFloat;
        #pragma unroll
        for (int i=0; i<numRegs; i++) penalty_here_array[i] = negInftyFloat;
        #pragma unroll //UNROLLHERE
        for (int i=0; i<numRegs; i++) F_here_array[i] = negInftyFloat;
        if (threadIdx.x % group_size == 0) {
            penalty_left = 0;
            penalty_diag = 0;
            #pragma unroll
            for (int i=0; i<numRegs; i++) penalty_here_array[i] = 0;
        }
        if (threadIdx.x % group_size == 1) {
            penalty_left = 0;
        }
    }

    __device__
    void init_local_score_profile_BLOSUM62(int offset_isc, int (&subject)[numRegs], 
        const char* const devS0, const int length_S0
    ) const{
        if (!offset_isc) {
            for (int i=threadIdx.x; i<21*21; i+=32) shared_BLOSUM62[i/21][i%21]=cBLOSUM62_dev[i];
            __syncwarp();
        }

        #pragma unroll //UNROLLHERE
        for (int i=0; i<numRegs; i++) {

            if (offset_isc+numRegs*(threadIdx.x%group_size)+i >= length_S0) subject[i] = 1; // 20;
            else{                
                subject[i] = devS0[offset_isc+numRegs*(threadIdx.x%group_size)+i];
            }
        }

    }

    __device__
    void shuffle_query(char new_letter, char& query_letter) const{
        query_letter = __shfl_up_sync(0xFFFFFFFF, query_letter, 1, 32);
        const int group_id = threadIdx.x % group_size;
        if (!group_id) query_letter = new_letter;
    }

    __device__
    void shuffle_new_query(char4& new_query_letter4) const{
        const int temp = __shfl_down_sync(0xFFFFFFFF, *((int*)(&new_query_letter4)), 1, 32);
        new_query_letter4 = *((char4*)(&temp));
    }

    __device__
    void shuffle_affine_penalty(
        float new_penalty_left, float new_E_left, float& E, 
        float penalty_here31, float& penalty_diag, float& penalty_left
    ) const{
        penalty_diag = penalty_left;
        penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here31, 1, 32);
        E = __shfl_up_sync(0xFFFFFFFF, E, 1, 32);
        const int group_id = threadIdx.x % group_size;
        if (!group_id) {
            penalty_left = new_penalty_left;
            E = new_E_left;
        }
    }

    __device__
    void shuffle_H_E_temp_out(short2& H_temp_out, short2& E_temp_out) const{
        const uint32_t temp = __shfl_down_sync(0xFFFFFFFF, *((int*)(&H_temp_out)), 1, 32);
        H_temp_out = *((short2*)(&temp));
        const uint32_t temp2 = __shfl_down_sync(0xFFFFFFFF, *((int*)(&E_temp_out)), 1, 32);
        E_temp_out = *((short2*)(&temp2));
    }

    __device__
    void shuffle_H_E_temp_in(short2& H_temp_in, short2& E_temp_in) const{
        const uint32_t temp = __shfl_down_sync(0xFFFFFFFF, *((int*)(&H_temp_in)), 1, 32);
        H_temp_in = *((short2*)(&temp));
        const uint32_t temp2 = __shfl_down_sync(0xFFFFFFFF, *((int*)(&E_temp_in)), 1, 32);
        E_temp_in = *((short2*)(&temp2));
    }

    __device__
    void set_H_E_temp_out_x(float penalty_here31, float E, short2& H_temp_out, short2& E_temp_out) const{
        if (threadIdx.x == 31) {
            H_temp_out.x = penalty_here31;
            E_temp_out.x = E;
        }
    };

    __device__
    void set_H_E_temp_out_y(float penalty_here31, float E, short2& H_temp_out, short2& E_temp_out) const{
        if (threadIdx.x == 31) {
            H_temp_out.y = penalty_here31;
            E_temp_out.y = E;
        }
    };



    __device__
    void computeFirstPass(float& maximum, const char* const devS0, const int length_S0
    ) const{
        // FIRST PASS (of many passes)
        // Note first pass has always full seqeunce length

        int counter = 1;
        char query_letter = 20;
        char4 new_query_letter4 = constantQuery4[threadIdx.x%group_size];
        if (threadIdx.x % group_size== 0) query_letter = new_query_letter4.x;


        const size_t base_3 = size_t(blockIdx.x)*size_t(length_2);
        short2* const devTempHcol = (&devTempHcol2[base_3]);
        short2* const devTempEcol = (&devTempEcol2[base_3]);

        const int group_id = threadIdx.x % group_size;
        int offset = group_id + group_size;
        int offset_out = group_id;


        float E = negInftyFloat;
        float penalty_here31;
        float penalty_diag;
        float penalty_left;
        short2 H_temp_out;
        short2 E_temp_out;
        int subject[numRegs];
        float penalty_here_array[numRegs];
        float F_here_array[numRegs];
       
        init_penalties_local(0, penalty_diag, penalty_left, penalty_here_array, F_here_array);
        init_local_score_profile_BLOSUM62(0, subject, devS0, length_S0);
        initial_calc32_local_affine_float(0, query_letter, E, penalty_here31, penalty_diag, penalty_left, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.y, query_letter);
        shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

        //shuffle_max();
        calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.z, query_letter);
        shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

        //shuffle_max();
        calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.w, query_letter);
        shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
        shuffle_new_query(new_query_letter4);
        counter++;

        for (int k = 4; k <= 28; k+=4) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.y, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.z, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.w, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_new_query(new_query_letter4);
            counter++;
        }

        for (int k = 32; k <= length_2+28; k+=4) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);            
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_H_E_temp_out(H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.y, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.z, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);

            if ((counter+8)%16 == 0 && counter > 8) {
                checkHEindex(offset_out, __LINE__);
                devTempHcol[offset_out]=H_temp_out;
                devTempEcol[offset_out]=E_temp_out;
                offset_out += group_size;
            }

            shuffle_H_E_temp_out(H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.w, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_new_query(new_query_letter4);
            if (counter%group_size == 0) {
                new_query_letter4 = constantQuery4[offset];
                offset += group_size;
            }
            counter++;
        }
        if (length_2 % 4 == 0) {
            const int temp1 = __shfl_up_sync(0xFFFFFFFF, *((int*)(&H_temp_out)), 1, 32);
            H_temp_out = *((short2*)(&temp1));
            const int temp2 = __shfl_up_sync(0xFFFFFFFF, *((int*)(&E_temp_out)), 1, 32);
            E_temp_out = *((short2*)(&temp2));
        }

        if (length_2%4 == 1) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
        }

        if (length_2%4 == 2) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
        }
        if (length_2%4 == 3) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_H_E_temp_out(H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.y, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
        }
        const int final_out = length_2 % 64;
        const int from_thread_id = 32 - ((final_out+1)/2);

        //printf("tid %d, offset_out %d, from_thread_id %d\n", threadIdx.x, offset_out, from_thread_id);
        if (threadIdx.x>=from_thread_id) {
            checkHEindex(offset_out-from_thread_id, __LINE__);
            devTempHcol[offset_out-from_thread_id]=H_temp_out;
            devTempEcol[offset_out-from_thread_id]=E_temp_out;
        }
    }

    __device__
    void computeMiddlePass(int pass, float& maximum, const char* const devS0, const int length_S0
    ) const{
        int counter = 1;
        char query_letter = 20;
        char4 new_query_letter4 = constantQuery4[threadIdx.x%group_size];
        if (threadIdx.x % group_size== 0) query_letter = new_query_letter4.x;


        const size_t base_3 = size_t(blockIdx.x)*size_t(length_2);
        short2* const devTempHcol = (&devTempHcol2[base_3]);
        short2* const devTempEcol = (&devTempEcol2[base_3]);

        const int group_id = threadIdx.x % group_size;
        int offset = group_id + group_size;
        int offset_out = group_id;
        int offset_in = group_id;
        checkHEindex(offset_in, __LINE__);
        short2 H_temp_in = devTempHcol[offset_in];
        short2 E_temp_in = devTempEcol[offset_in];
        offset_in += group_size;

        float E = negInftyFloat;
        float penalty_here31;
        float penalty_diag;
        float penalty_left;
        short2 H_temp_out;
        short2 E_temp_out;
        int subject[numRegs];
        float penalty_here_array[numRegs];
        float F_here_array[numRegs];

        init_penalties_local(gap_open+(pass*32*numRegs-1)*gap_extend, penalty_diag, penalty_left, penalty_here_array, F_here_array);
        init_local_score_profile_BLOSUM62(pass*(32*numRegs), subject, devS0, length_S0);

        if (!group_id) {
            penalty_left = H_temp_in.x;
            E = E_temp_in.x;
        }


        initial_calc32_local_affine_float(1, query_letter, E, penalty_here31, penalty_diag, penalty_left, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.y, query_letter);
        shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
        shuffle_H_E_temp_in(H_temp_in, E_temp_in);

        //shuffle_max();
        calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.z, query_letter);
        shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

        //shuffle_max();
        calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.w, query_letter);
        shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
        shuffle_H_E_temp_in(H_temp_in, E_temp_in);
        shuffle_new_query(new_query_letter4);
        counter++;

        for (int k = 4; k <= 28; k+=4) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.y, query_letter);
            shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_H_E_temp_in(H_temp_in, E_temp_in);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.z, query_letter);
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.w, query_letter);
            shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_new_query(new_query_letter4);
            shuffle_H_E_temp_in(H_temp_in, E_temp_in);

            counter++;
        }
        for (int k = 32; k <= length_2+28; k+=4) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_H_E_temp_out(H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.y, query_letter);
            shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_H_E_temp_in(H_temp_in, E_temp_in);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.z, query_letter);
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);

            if ((counter+8)%16 == 0 && counter > 8) {
                checkHEindex(offset_out, __LINE__);
                devTempHcol[offset_out]=H_temp_out;
                devTempEcol[offset_out]=E_temp_out;
                offset_out += group_size;
            }
            shuffle_H_E_temp_out(H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.w, query_letter);
            shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_new_query(new_query_letter4);
            if (counter%group_size == 0) {
                new_query_letter4 = constantQuery4[offset];
                offset += group_size;
            }
            shuffle_H_E_temp_in(H_temp_in, E_temp_in);
            if (counter%16 == 0) {
                checkHEindex(offset_in, __LINE__);
                H_temp_in = devTempHcol[offset_in];
                E_temp_in = devTempEcol[offset_in];
                offset_in += group_size;
            }
            counter++;
        }

        if (length_2 % 4 == 0) {
            const int temp1 = __shfl_up_sync(0xFFFFFFFF, *((int*)(&H_temp_out)), 1, 32);
            H_temp_out = *((short2*)(&temp1));
            const int temp2 = __shfl_up_sync(0xFFFFFFFF, *((int*)(&E_temp_out)), 1, 32);
            E_temp_out = *((short2*)(&temp2));
        }
        if (length_2 % 4 == 1) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            set_H_E_temp_out_y(penalty_here31, E, H_temp_out, E_temp_out);
        }        

        if (length_2%4 == 2) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(penalty_here31, E, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(E, penalty_here31, H_temp_out, E_temp_out);
        }
        if (length_2%4 == 3) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(E, penalty_here31, H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.x, query_letter);
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_y(E, penalty_here31, H_temp_out, E_temp_out);
            shuffle_H_E_temp_out(H_temp_out, E_temp_out);
            shuffle_query(new_query_letter4.y, query_letter);
            shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);

            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            set_H_E_temp_out_x(E, penalty_here31, H_temp_out, E_temp_out);
            set_H_E_temp_out_y(E, penalty_here31, H_temp_out, E_temp_out);
        }
        const int final_out = length_2 % 64;
        const int from_thread_id = 32 - ((final_out+1)/2);

        if (threadIdx.x>=from_thread_id) {
            checkHEindex(offset_out-from_thread_id, __LINE__);
            devTempHcol[offset_out-from_thread_id]=H_temp_out;
            devTempEcol[offset_out-from_thread_id]=E_temp_out;
        }
    }

    __device__ 
    void computeFinalPass(int passes, float& maximum, const char* const devS0, const int length_S0
    ) const{
        int counter = 1;
        char query_letter = 20;
        char4 new_query_letter4 = constantQuery4[threadIdx.x%group_size];
        if (threadIdx.x % group_size== 0) query_letter = new_query_letter4.x;


        const size_t base_3 = size_t(blockIdx.x)*size_t(length_2);
        short2* const devTempHcol = (&devTempHcol2[base_3]);
        short2* const devTempEcol = (&devTempEcol2[base_3]);

        const int group_id = threadIdx.x % group_size;
        int offset = group_id + group_size;
        int offset_in = group_id;
        checkHEindex(offset_in, __LINE__);
        short2 H_temp_in = devTempHcol[offset_in];
        short2 E_temp_in = devTempEcol[offset_in];
        offset_in += group_size;

        const uint32_t thread_result = ((length_S0-1)%(32*numRegs))/numRegs;

        float E = negInftyFloat;
        float penalty_here31;
        float penalty_diag;
        float penalty_left;
        int subject[numRegs];
        float penalty_here_array[numRegs];
        float F_here_array[numRegs];

        init_penalties_local(gap_open+((passes-1)*32*numRegs-1)*gap_extend, penalty_diag, penalty_left, penalty_here_array, F_here_array);
        init_local_score_profile_BLOSUM62((passes-1)*(32*numRegs), subject, devS0, length_S0);
        //copy_H_E_temp_in();
        if (!group_id) {
            penalty_left = H_temp_in.x;
            E = E_temp_in.x;
        }

        initial_calc32_local_affine_float(1, query_letter, E, penalty_here31, penalty_diag, penalty_left, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.y, query_letter);
        shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
        shuffle_H_E_temp_in(H_temp_in, E_temp_in);
        if (length_2+thread_result >=2) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.z, query_letter);
            //copy_H_E_temp_in();
            shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);
        }

        if (length_2+thread_result >=3) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.w, query_letter);
            //copy_H_E_temp_in();
            shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_H_E_temp_in(H_temp_in, E_temp_in);
            shuffle_new_query(new_query_letter4);
            counter++;
        }
        if (length_2+thread_result >=4) {
            int k;
            //for (k = 5; k < lane_2+thread_result-2; k+=4) {
            for (k = 4; k <= length_2+(thread_result-3); k+=4) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);

                shuffle_query(new_query_letter4.x, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.y, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
                shuffle_H_E_temp_in(H_temp_in, E_temp_in);

                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.z, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);

                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.w, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
                shuffle_new_query(new_query_letter4);
                if (counter%group_size == 0) {
                    new_query_letter4 = constantQuery4[offset];
                    offset += group_size;
                }
                shuffle_H_E_temp_in(H_temp_in, E_temp_in);
                if (counter%16 == 0) {
                    checkHEindex(offset_in, __LINE__);
                    H_temp_in = devTempHcol[offset_in];
                    E_temp_in = devTempEcol[offset_in];
                    offset_in += group_size;
                }
                counter++;
            }

            if ((k-1)-(length_2+thread_result) > 0) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.x, query_letter);
                shuffle_affine_penalty(H_temp_in.x, E_temp_in.x, E, penalty_here31, penalty_diag, penalty_left);
                k++;
            }


            if ((k-1)-(length_2+thread_result) > 0) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.y, query_letter);
                shuffle_affine_penalty(H_temp_in.y, E_temp_in.y, E, penalty_here31, penalty_diag, penalty_left);
                shuffle_H_E_temp_in(H_temp_in, E_temp_in);
                k++;
            }

            if ((k-1)-(length_2+thread_result) > 0) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            }
        }
    }

    __device__ 
    void computeSinglePass(float& maximum, const char* const devS0, const int length_S0
    ) const{
        int counter = 1;
        char query_letter = 20;
        char4 new_query_letter4 = constantQuery4[threadIdx.x%group_size];
        if (threadIdx.x % group_size== 0) query_letter = new_query_letter4.x;

        const int group_id = threadIdx.x % group_size;
        int offset = group_id + group_size;
        int offset_in = group_id;
        checkHEindex(offset_in, __LINE__);
        offset_in += group_size;

        const uint32_t thread_result = ((length_S0-1)%(32*numRegs))/numRegs;

        float E = negInftyFloat;
        float penalty_here31;
        float penalty_diag;
        float penalty_left;
        int subject[numRegs];
        float penalty_here_array[numRegs];
        float F_here_array[numRegs];

        init_penalties_local(0, penalty_diag, penalty_left, penalty_here_array, F_here_array);
        init_local_score_profile_BLOSUM62(0, subject, devS0, length_S0);

        initial_calc32_local_affine_float(0, query_letter, E, penalty_here31, penalty_diag, penalty_left, maximum, subject, penalty_here_array, F_here_array);
        shuffle_query(new_query_letter4.y, query_letter);
        shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

        if (length_2+thread_result >=2) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.z, query_letter);
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
        }

        if (length_2+thread_result >=3) {
            //shuffle_max();
            calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            shuffle_query(new_query_letter4.w, query_letter);
            //copy_H_E_temp_in();
            shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
            shuffle_new_query(new_query_letter4);
            counter++;
        }
        if (length_2+thread_result >=4) {
            int k;
            //for (k = 5; k < lane_2+thread_result-2; k+=4) {
            for (k = 4; k <= length_2+(thread_result-3); k+=4) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);

                shuffle_query(new_query_letter4.x, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.y, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.z, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);

                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.w, query_letter);
                //copy_H_E_temp_in();
                shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
                shuffle_new_query(new_query_letter4);
                if (counter%group_size == 0) {
                    new_query_letter4 = constantQuery4[offset];
                    offset += group_size;
                }
                counter++;
            }

            if ((k-1)-(length_2+thread_result) > 0) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.x, query_letter);
                shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
                k++;
            }


            if ((k-1)-(length_2+thread_result) > 0) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
                shuffle_query(new_query_letter4.y, query_letter);
                shuffle_affine_penalty(0.f, negInftyFloat, E, penalty_here31, penalty_diag, penalty_left);
                k++;
            }

            if ((k-1)-(length_2+thread_result) > 0) {
                //shuffle_max();
                calc32_local_affine_float(query_letter, E, penalty_here31, penalty_diag, maximum, subject, penalty_here_array, F_here_array);
            }
        }
    }

    template<class ScoreOutputIterator>
    __device__
    void compute(
        ScoreOutputIterator const devAlignmentScores
    ) const{


        const int length_S0 = devLengths[d_positions_of_selected_lengths[blockIdx.x]];
        const size_t base_S0 = devOffsets[d_positions_of_selected_lengths[blockIdx.x]]-devOffsets[0];

        const char* const devS0 = &devChars[base_S0];

        const int passes = (length_S0 + (group_size*numRegs) - 1) / (group_size*numRegs);

        float maximum = 0.f;

        if(passes == 1){
            computeSinglePass(maximum, devS0, length_S0);
        }else{

            computeFirstPass(maximum, devS0, length_S0);

            for (int pass = 1; pass < passes-1; pass++) {
                computeMiddlePass(pass, maximum, devS0, length_S0);
            }

            computeFinalPass(passes, maximum, devS0, length_S0);
        }

        for (int offset=group_size/2; offset>0; offset/=2){
            maximum = max(maximum,__shfl_down_sync(0xFFFFFFFF,maximum,offset,group_size));
        }

        const int group_id = threadIdx.x % group_size;
        if (!group_id) {
            devAlignmentScores[d_positions_of_selected_lengths[blockIdx.x]] = maximum;
        }
    }
};



// Needleman-Wunsch (NW): global alignment with linear gap penalty
// numRegs values per thread
// uses a single warp per CUDA thread block;
// every groupsize threads computes an alignmen score
template <int numRegs, class ScoreOutputIterator, class PositionsIterator> 
__launch_bounds__(32,16)
//__launch_bounds__(32)
__global__
void NW_local_affine_read4_float_query_Protein_new(
    __grid_constant__ const char * const devChars,
    __grid_constant__ ScoreOutputIterator const devAlignmentScores,
    __grid_constant__ short2 * const devTempHcol2,
    __grid_constant__ short2 * const devTempEcol2,
    __grid_constant__ const size_t* const devOffsets,
    __grid_constant__ const size_t* const devLengths,
    __grid_constant__ PositionsIterator const d_positions_of_selected_lengths,
    __grid_constant__ const int length_2,
    __grid_constant__ const float gap_open,
    __grid_constant__ const float gap_extend
) {
    using Processor = ManyPassFloat<numRegs, PositionsIterator>;

    __shared__ typename Processor::BLOSUM62_SMEM shared_BLOSUM62;

    Processor processor(
        shared_BLOSUM62,
        devChars,
        devTempHcol2,
        devTempEcol2,
        devOffsets,
        devLengths,
        d_positions_of_selected_lengths,
        length_2,
        gap_open,
        gap_extend
    );

    processor.compute(devAlignmentScores);
}


template <int numRegs, class ScoreOutputIterator> 
__launch_bounds__(1,1)
__global__
void launch_process_overflow_alignments_kernel_NW_local_affine_read4_float_query_Protein_new(
    __grid_constant__ const int* const d_overflow_number,
    __grid_constant__ short2* const d_temp,
    __grid_constant__ const size_t maxTempBytes,
    __grid_constant__ const char * const devChars,
    __grid_constant__ ScoreOutputIterator const devAlignmentScores,
    __grid_constant__ const size_t* const devOffsets,
    __grid_constant__ const size_t* const devLengths,
    __grid_constant__ const size_t* const d_positions_of_selected_lengths,
    __grid_constant__ const int queryLength,
    __grid_constant__ const float gap_open,
    __grid_constant__ const float gap_extend
){
    const int numOverflow = *d_overflow_number;
    if(numOverflow > 0){
        // printf("numOverflow %d\n", numOverflow);
        // for(int i = 0; i < numOverflow; i++){
        //     printf("%lu ", d_positions_of_selected_lengths[i]);

        //     size_t p = d_positions_of_selected_lengths[i];
        //     size_t offset = devOffsets[p];
        //     int length = devLengths[p];
        //     printf("length %d\n", length);
        //     for(int k = 0; k < length; k++){
        //         printf("%d ", int(devChars[offset + k]));
        //     }
        //     printf("\n");
        // }
        // printf("\n");
        // printf("next 5");
        // for(int i = numOverflow; i < numOverflow + 5; i++){
        //     printf("%lu ", d_positions_of_selected_lengths[i]);
        // }
        // printf("\n");
    
        const size_t tempBytesPerSubjectPerBuffer = sizeof(short2) * queryLength;
        const size_t maxSubjectsPerIteration = std::min(size_t(numOverflow), maxTempBytes / (tempBytesPerSubjectPerBuffer * 2));

        short2* d_tempHcol2 = d_temp;
        short2* d_tempEcol2 = (short2*)(((char*)d_tempHcol2) + maxSubjectsPerIteration * tempBytesPerSubjectPerBuffer);

        const int numIters =  SDIV(numOverflow, maxSubjectsPerIteration);
        for(int iter = 0; iter < numIters; iter++){
            const size_t begin = iter * maxSubjectsPerIteration;
            const size_t end = iter < numIters-1 ? (iter+1) * maxSubjectsPerIteration : numOverflow;
            const size_t num = end - begin;

            cudaMemsetAsync(d_temp, 0, tempBytesPerSubjectPerBuffer * 2 * num, 0);

            // cudaMemsetAsync(d_tempHcol2, 0, tempBytesPerSubjectPerBuffer * num, 0);
            // cudaMemsetAsync(d_tempEcol2, 0, tempBytesPerSubjectPerBuffer * num, 0);

            NW_local_affine_read4_float_query_Protein_new<12><<<num, 32>>>(
                devChars, 
                devAlignmentScores,
                d_tempHcol2, 
                d_tempEcol2, 
                devOffsets, 
                devLengths, 
                d_positions_of_selected_lengths + begin, 
                queryLength, 
                gap_open, 
                gap_extend
            );
        }
    }
}

#endif