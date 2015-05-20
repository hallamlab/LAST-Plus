//
// Created by david on 20/05/15.
//

#include "inputOutput.h"

std::queue<std::string> output;

void* writer_func(void* args) {
/*
	while(1) {
		sem_wait(sema_w);

		QUEUE* temp;

		sem_wait(sema_Q);
		cutnpaste_q(&temp, DONE_Q);
		sem_post(sema_Q);

		while(temp) {
			sem_wait(sema_R);

			thread_data* td = temp->td;
			unsigned int buffer = temp->buffer;
			num_writes++;

			for(j = 0; j < td->output_num_sequences[buffer]; j++) {
				writer_counter++;
				char *ptrc;

				if(td->aa_buffer[buffer][j][0]!=0) {
					ptrc=td->aa_buffer[buffer][j];

					while(*ptrc!='\0'){
						if(*ptrc=='\t') *ptrc='>';
						ptrc++;
					}
					fprintf(aa_outfile_fp, ">%s", td->aa_buffer[buffer][j]);
				}
				memset(td->aa_buffer[buffer][j], 0, STRINGLEN);
			}

			if (verbose) printf("INFO: Wrote results for thread %d, buffer %d.\n", td->id, buffer );

			if (output_meta) fclose(outfile_fp);
			if (output_dna) fclose(dna_outfile_fp);

			sem_post(sema_R);
			sem_post(td->sema_w);

			temp = temp->next;
		}

		if(num_reads_flag == 1 && writer_counter ==  read_counter)   {
			sem_post(stop_sema);
			break;
		}
	}
	fclose(aa_outfile_fp);
 */
}

void conductWork(){
/*
	while (stopped_at_fpos!=0) {
		sem_wait(sema_r);

		sem_wait(sema_Q);
		QUEUE* temp;
		cutnpaste_q(&temp, EMPTY_Q);
		sem_post(sema_Q);

		while(temp) {
			sem_wait(sema_R);

			stopped_at_fpos = read_seq_into_buffer(fp,  temp->td, temp->buffer);

			sem_post(sema_R);

			sem_post(temp->td->sema_r);
			temp = temp->next;
		}
	}
	CloseFASTA(fp);

	if (verbose) printf("INFO : Finished handing out all the work...\n");
	num_reads_flag =1;

	sem_wait(stop_sema);
 */
}


void* thread_func(void *_thread_datas) {
/*
	thread_data *td = (thread_data*)_thread_datas;
	unsigned int b = 0;
	unsigned int i;

	while(1) {
		sem_wait(td->sema_r);
		sem_wait(td->sema_w);

		sem_wait(counter_sema);
		viterbi_counter +=  td->input_num_sequences[b];
		sem_post(counter_sema);

		for (i=0; i < td->input_num_sequences[b]; i++) {
			unsigned int stringlength = strlen(td->input_buffer[b][i]);
			get_prob_from_cg(td->hmm, &train, td->input_buffer[b][i], stringlength);

			if(td->input_buffer[b][i] == 0 || td->input_head_buffer[b][i] == 0 ) {
				printf("%s\n",td->input_buffer[b][i]);
				printf("%s\n",td->input_head_buffer[b][i]);
			}

			if (td->input_buffer[b][i] != 0 && td->input_head_buffer[b][i] != 0 ) {
				memset(td->aa_buffer[b][i], 0, STRINGLEN );

				viterbi(td->hmm, td->input_buffer[b][i], td->output_buffer[b][i], td->aa_buffer[b][i], td->dna_buffer[b][i],
				        td->input_head_buffer[b][i], td->wholegenome, td->format, stringlength,
				        td->dna, td->dna1, td->dna_f, td->dna_f1, td->protein,
				        td->insert, td->c_delete, td->temp_str);

				td->acceptable_buffer[b][i] = 1;

				sem_wait(work_sema);
				work_counter++;
				sem_post(work_sema);
			}
		}
		td->output_num_sequences[b] = td->input_num_sequences[b];

		if (verbose) printf("INFO: Thread %d buffer %d done work on %d sequences!\n", td->id, b, td->input_num_sequences[b]);
		sem_wait(sema_Q);
		enqueue(td, b, EMPTY_Q);
		enqueue(td, b, DONE_Q);

		sem_post(sema_Q);
		sem_post(sema_r);
		sem_post(sema_w);

		b = (b + 1) % 2;
	}
	return (void*) 0;
 */
}
