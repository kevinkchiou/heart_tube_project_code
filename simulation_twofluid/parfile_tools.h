int trim_whitespace(char *out,char *in);
int locate_whitespace(char *str,int *idx);
char **split_string_whitespace(int *size,char *str);
void free_split_string(char **str,const unsigned int size);
