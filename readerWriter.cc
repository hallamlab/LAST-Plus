int read_seq_into_buffer(FASTAFILE* ffp, threadData &thread_data, unsigned int buf) {

  char *seq;
  char *name;
  int L;
  int status ;
  int count=0;

  while ( (count < MAX_SEQS_PER_BUFFER) && (status = ReadFASTA(ffp, &seq, &name, &L)) ==1 ) {
    //strcpy(thread_data->input_head_buffer[buf][count], name);
    //strcpy(thread_data->input_buffer[buf][count], seq);

    read_counter++;
    count++;
  }

  thread_data->input_num_sequences[buf] = count;
  read_counter1 += thread_data->input_num_sequences[buf];

  return count;

}

void* writer_func(void* args) {

  int j;
  //FILE* aa_outfile_fp = fopen(aa_file, "a");

  if(!aa_outfile_fp) {
    printf("ERROR: Could not open aa output file %s for writing!\n", aa_file);
    exit(0);
  }

  while(1) {
    SEM_WAIT(sema_w);
    QUEUE* temp;
    SEM_WAIT(sema_Q);
    cutnpaste_q(&temp, DONE_Q);
    SEM_POST(sema_Q);

    while(temp) {
      SEM_WAIT(sema_R);
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
      SEM_POST(sema_R);
      SEM_POST(td->sema_w);
      temp = temp->next;
    }

    if(num_reads_flag == 1 && writer_counter ==  read_counter)   {
      printf("TERMINATING\n");
      SEM_POST(stop_sema);
      break;
    }
  }
  fclose(aa_outfile_fp);
}

