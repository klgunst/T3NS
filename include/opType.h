#ifndef OPTYPE_H
# define OPTYPE_H

struct opType {
        int *begin_opType; /* ordered like 0,1,2,3,4 
                            * and normal/complimentary */
        int ***tags_opType;
};
int amount_opType(const struct opType * const ops, const int nrop, 
                  const char t);

int range_opType(int * const begin, int * const end, 
                  const struct opType * const ops, const int nrop);

int id_opType(const struct opType * const ops, const char c);

void get_opType(struct opType * const ops, const int bond, const int is_left);

void get_opType_site(struct opType * const ops, const int psite);

void get_unity_opType(struct opType * const ops);

double interactval(const int ids[3], const struct opType ops[3], const char t);

void init_opType_array(void);

void destroy_opType(struct opType * const ops, const int bond, 
                    const int is_left);

void symsec_of_operators(int ** const list_of_ss, const int bond, 
                                const int is_left);

void get_todo_and_creator_array(const int (**todo)[2], 
                                const int (**creator_array)[2][3][2]);

int get_combine_array(const int (**array)[3]);

int fillin_nr_basetags(int nr_basetags[5][2]);

int loop_dof(const int nr, const int creator[nr], const int position[nr], 
             int dof[nr], const char t, const int bond, const int is_left);

void opType_get_string_of_rops(char buffer[], const int ropsindex, 
                               const int bond, const int is_left);

void opType_get_string_of_siteops(char buffer[], const int siteid, 
                                  const int site);

int opType_symsec_siteop(const int siteoperator, const int site);
#endif
