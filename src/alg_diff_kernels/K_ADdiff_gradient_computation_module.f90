module ADdiff_gradient_computation_module
implicit none
contains


!  Differentiation of compute_grad_scalarfield_gg in forward (tangent) mode:
!   variations   of useful results: field_grad
!   with respect to varying inputs: field
!   RW status of diff variables: field:in field_grad:out
!=======================================================================
!============ GREEN GAUSS GRADIENT COMPUTATION SUBROUTINE ==============
!=======================================================================
	!=== FIELD SUBS: ===
  PURE SUBROUTINE COMPUTE_GRAD_SCALARFIELD_GG_D(fieldd, field_gradd, &
&   face_left_cell, face_right_cell, face_normal, face_area&
&   , cell_volume, face_weight, ncells, nfaces, nbfaces, dim, &
&   use_ghost_cells)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncells, nfaces, nbfaces, dim
    LOGICAL, INTENT(IN) :: use_ghost_cells
    INTEGER, INTENT(IN), CONTIGUOUS :: face_left_cell(:), &
&   face_right_cell(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: face_normal(:, :), face_area(:), &
&   cell_volume(:), face_weight(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: fieldd(:)
    REAL*8, INTENT(INOUT), CONTIGUOUS :: field_gradd(:, :)
    INTEGER :: n_loop_faces
    INTEGER :: i, left_cell, right_cell
    REAL*8 :: w, inv_vol
    REAL*8 :: avg_face_valued, vald(dim)
	
	field_gradd = 0.0_8
    n_loop_faces = nfaces - nbfaces
    IF (use_ghost_cells) THEN
      n_loop_faces = nfaces
    END IF
!=== INTERNAL FACES: ===
    DO i=1,n_loop_faces
      left_cell = face_left_cell(i)
      right_cell = face_right_cell(i)
      w = face_weight(i)
      avg_face_valued = w*fieldd(right_cell) + (1.d0-w)*fieldd(left_cell&
&       )
      vald = face_normal(:, i)*face_area(i)*avg_face_valued
      field_gradd(:, left_cell) = field_gradd(:, left_cell) + vald
      field_gradd(:, right_cell) = field_gradd(:, right_cell) - vald
    END DO
    IF (.NOT.use_ghost_cells) THEN
!=== UPDATING ACCORDING BOUNDARY FACES VALUES: ===
      DO i=nfaces-nbfaces+1,nfaces
        left_cell = face_left_cell(i)
        right_cell = face_right_cell(i)
        avg_face_valued = fieldd(right_cell)
        vald = face_normal(:, i)*face_area(i)*avg_face_valued
        field_gradd(:, left_cell) = field_gradd(:, left_cell) + vald
      END DO
    END IF
!=== DEVIDING BY CELL VOLUME: ===
    DO i=1,ncells
      inv_vol = 1.d0/cell_volume(i)
      field_gradd(:, i) = inv_vol*field_gradd(:, i)
    END DO
  END SUBROUTINE

!  Differentiation of compute_grad_vectorfield_gg in forward (tangent) mode:
!   variations   of useful results: field_grad
!   with respect to varying inputs: field
!   RW status of diff variables: field:in field_grad:out
  PURE SUBROUTINE COMPUTE_GRAD_VECTORFIELD_GG_D(fieldd, field_gradd&
&   , face_left_cell, face_right_cell, face_normal, face_area&
&   , cell_volume, face_weight, ncells, nfaces, nbfaces, dim, nvars, &
&   use_ghost_cells)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncells, nfaces, nbfaces, dim, nvars
    LOGICAL, INTENT(IN) :: use_ghost_cells
    INTEGER, INTENT(IN), CONTIGUOUS :: face_left_cell(:), &
&   face_right_cell(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: face_normal(:, :), face_area(:), &
&   cell_volume(:), face_weight(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: fieldd(:, :)
    REAL*8, INTENT(INOUT), CONTIGUOUS :: field_gradd(:, :)
    INTEGER :: n_loop_faces
    INTEGER :: i, j, d, offset
    INTEGER :: left_cell, right_cell
    REAL*8 :: w, inv_vol
    REAL*8 :: avg_face_valued(nvars), vald
    field_gradd = 0.0_8
    n_loop_faces = nfaces - nbfaces
    IF (use_ghost_cells) THEN
      n_loop_faces = nfaces
    END IF
!=== INTERNAL FACES: ===
    DO i=1,n_loop_faces
      left_cell = face_left_cell(i)
      right_cell = face_right_cell(i)
      w = face_weight(i)
      avg_face_valued = w*fieldd(:, right_cell) + (1.d0-w)*fieldd(:, &
&       left_cell)
      DO j=1,nvars
        offset = (j-1)*dim
        DO d=1,dim
          vald = face_normal(d, i)*face_area(i)*avg_face_valued(j)
          field_gradd(d+offset, left_cell) = field_gradd(d+offset, &
&           left_cell) + vald
          field_gradd(d+offset, right_cell) = field_gradd(d+offset, &
&           right_cell) - vald
        END DO
      END DO
    END DO
    IF (.NOT.use_ghost_cells) THEN
!=== UPDATING ACCORDING BOUNDARY FACES VALUES: ===
      DO i=nfaces-nbfaces+1,nfaces
        left_cell = face_left_cell(i)
        right_cell = face_right_cell(i)
        avg_face_valued = fieldd(:, right_cell)
        DO j=1,nvars
          offset = (j-1)*dim
          DO d=1,dim
            vald = face_normal(d, i)*face_area(i)*avg_face_valued(j)
            field_gradd(d+offset, left_cell) = field_gradd(d+offset, &
&             left_cell) + vald
          END DO
        END DO
      END DO
    END IF
!=== DEVIDING BY CELL VOLUME: ===
    DO i=1,ncells
      inv_vol = 1.d0/cell_volume(i)
      field_gradd(:, i) = inv_vol*field_gradd(:, i)
    END DO
  END SUBROUTINE


!=======================================================================
!============ LEAST SQUARES GRADIENT COMPUTATION SUBROUTINE ============
!=======================================================================
	!=== FIELD SUBS: ===

!  Differentiation of compute_grad_scalarfield_lsq in forward (tangent) mode:
!   variations   of useful results: field_grad
!   with respect to varying inputs: field
!   RW status of diff variables: field:in field_grad:out
  PURE SUBROUTINE COMPUTE_GRAD_SCALARFIELD_LSQ_D(fieldd, field_gradd&
&   , ncells, dim, cell_faces_ptr, cell_faces, face_left_cell&
&   , face_right_cell, lsq_w)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncells, dim
    INTEGER, INTENT(IN), CONTIGUOUS :: cell_faces_ptr(:), cell_faces(:)&
&   , face_left_cell(:), face_right_cell(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: fieldd(:) 						
    REAL*8, INTENT(IN), CONTIGUOUS :: lsq_w(:, :)						! (dim, size(cell_faces))
    REAL*8, INTENT(INOUT), CONTIGUOUS :: field_gradd(:, :)
    INTEGER :: cell_idx, f, face_idx, neighbor, pos1, pos2
    REAL*8 :: weight(3)
    REAL*8 :: dphid
    field_gradd = 0.0_8
    DO cell_idx=1,ncells
      pos1 = cell_faces_ptr(cell_idx)
      pos2 = cell_faces_ptr(cell_idx+1) - 1
      DO f=pos1,pos2
        face_idx = cell_faces(f)
        IF (face_left_cell(face_idx) .EQ. cell_idx) THEN
          neighbor = face_right_cell(face_idx)
        ELSE
          neighbor = face_left_cell(face_idx)
        END IF
        weight(1:dim) = lsq_w(1:dim, f)
        dphid = fieldd(neighbor) - fieldd(cell_idx)
        field_gradd(:, cell_idx) = field_gradd(:, cell_idx) + weight(1:&
&         dim)*dphid
      END DO
    END DO
  END SUBROUTINE


!  Differentiation of compute_grad_vectorfield_lsq in forward (tangent) mode:
!   variations   of useful results: field_grad
!   with respect to varying inputs: field
!   RW status of diff variables: field:in field_grad:out
  PURE SUBROUTINE COMPUTE_GRAD_VECTORFIELD_LSQ_D(fieldd, field_gradd&
&   , ncells, nvars, dim, cell_faces_ptr, cell_faces, &
&   face_left_cell, face_right_cell, lsq_w)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ncells, nvars, dim
    INTEGER, INTENT(IN), CONTIGUOUS :: cell_faces_ptr(:), cell_faces(:)&
&   , face_left_cell(:), face_right_cell(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: fieldd(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: lsq_w(:, :)						! (dim, size(cell_faces))
    REAL*8, INTENT(INOUT), CONTIGUOUS :: field_gradd(:, :)
    INTEGER :: cell_idx, f, face_idx, neighbor, v, d, pos1, pos2, offset
    REAL*8 ::  weight(3)
    REAL*8 :: dphid
    field_gradd = 0.0_8
    DO cell_idx=1,ncells
      pos1 = cell_faces_ptr(cell_idx)
      pos2 = cell_faces_ptr(cell_idx+1) - 1
      DO f=pos1,pos2
        face_idx = cell_faces(f)
        IF (face_left_cell(face_idx) .EQ. cell_idx) THEN
          neighbor = face_right_cell(face_idx)
        ELSE
          neighbor = face_left_cell(face_idx)
        END IF
        weight(1:dim) = lsq_w(1:dim, f)
        DO v=1,nvars
          dphid = fieldd(v, neighbor) - fieldd(v, cell_idx)
          offset = (v-1)*dim
          DO d=1,dim
            field_gradd(d+offset, cell_idx) = field_gradd(d+offset, &
&             cell_idx) + weight(d)*dphid
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE

end module
