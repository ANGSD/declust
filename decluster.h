/*
 *  * Macro:[ASSERT]
 *   * shortcut to evaluate an expression, works the same way as the C-macro assert
 *    */
#define ASSERT(expr) if (!(expr)) {fprintf(stderr,"\n\n*******\n[ERROR](%s:%d) %s\n*******\n",__FILE__,__LINE__,#expr);exit(1);}

#define DEBUG(expr) if (DEBUG_MODE==1) {fprintf(stderr,"HEY");}
// {expr;}
