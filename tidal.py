# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 17:25:06 2016

@author: bperfect
"""

#Western edge
fac=TANH((tdays(ng)-dstart)/1.0_r8)
omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
minor=0.0143_r8+(0.0143_r8+0.010_r8)/REAL(Iend+1,r8)
major=0.1144_r8+(0.1144_r8-0.013_r8)/REAL(Iend+1,r8)
phase=(318.0_r8+(318.0_r8-355.0_r8)/REAL(Iend+1,r8))*deg2rad
angle=(125.0_r8+(125.0_r8- 25.0_r8)/REAL(Iend+1,r8))*deg2rad
DO j=JstrT,JendT
    val=0.5_r8*(angler(Istr-1,j)+angler(Istr,j))
    
ubar_west(j)=0.5_r8+fac*(major*COS(angle-val)*COS(omega-phase) - minor*SIN(angle-val)*SIN(omega-phase))
        END DO
        DO j=JstrP,JendT
          val=0.5_r8*(angler(Istr-1,j-1)+angler(Istr-1,j))
          BOUNDARY(ng)%vbar_west(j)=fac*(major*SIN(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         COS(omega-phase))
        END DO
      END IF

#Iend is supposed to be the number of points in the I direction
#angler is an input matrix that might be the angle doodles that only matter in a global model

#Eastern edge
#fac=TANH((tdays(ng)-dstart)/1.0_r8)
#omega=2.0_r8*pi*time(ng)/(12.42_r8*3600.0_r8)  !  M2 Tide period
        minor=0.0143_r8+(0.0143_r8+0.010_r8)
        major=0.1144_r8+(0.1144_r8-0.013_r8)
        phase=(318.0_r8+(318.0_r8-355.0_r8))*deg2rad
        angle=(125.0_r8+(125.0_r8- 25.0_r8))*deg2rad
        DO j=JstrT,JendT
          val=0.5_r8*(angler(Iend,j)+angler(Iend+1,j))
          BOUNDARY(ng)%ubar_east(j)=0.5_r8+fac*(major*COS(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         SIN(omega-phase))
        END DO
        DO j=JstrP,JendT
          val=0.5_r8*(angler(Iend+1,j-1)+angler(Iend+1,j))
          BOUNDARY(ng)%vbar_east(j)=fac*(major*SIN(angle-val)*          &
     &                                         COS(omega-phase)-        &
     &                                   minor*SIN(angle-val)*          &
     &                                         COS(omega-phase))
        END DO
      END IF
