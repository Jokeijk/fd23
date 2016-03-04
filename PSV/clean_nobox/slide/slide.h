#ifndef _SLIDE_H_
#define _SLIDE_H_
class PARAM;
class SLIDE{
  public:
  int current_step;
  int nxblock;
  int nzblock;
  int* first_arrive_step;
  int* last_step;
  SLIDE(PARAM & param);
  int usepunch;
};
#endif
