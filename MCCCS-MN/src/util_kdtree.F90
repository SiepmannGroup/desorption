MODULE util_kdtree
    use var_type, only:dp
    use util_runtime,only:err_exit
    use util_random,only:random
    use sim_system,only:interval,tree_node,tree
    implicit none
    private

    ! public subroutines
    public :: init_tree, insert_node, range_search, find_min, find_max, delete_node, &
             check_tree_coord, check_tree_cube, search_node_in_tree, update_cube, &
             init_tree_coordinates, update_tree_height, vector_dist_square, &
             empty_tree, allocate_kdtree, construct_kdtree, update_box_kdtree, read_kdtree, &
             allocate_nxyz, scale_kdtree

    ! global variables
    real, public :: rmin_square, rcut_square
    real :: coord_replaced(3) !< used in the deletion

contains

    ! comparison in kd-tree structure
    ! a_i, a_j: the first dimension of two coordinates to be compared
    ! b_i, b_j: the second ...
    ! c_i, c_j: the third ...
    ! return: whether the coordinate of bead i is "smaller" than that of bead j
    ! compare a dimension first, if ties, compare b dimension, if ties, compare c dimension
    function smallerThan(a_i, a_j, b_i, b_j, c_i, c_j)
        real :: a_i, a_j, b_i, b_j, c_i, c_j
        logical :: smallerThan

        if ((a_i - a_j) .lt. -1e-6) then
            smallerThan = .true.
        else if ((a_i - a_j) .gt. 1e-6) then
            smallerThan = .false.
        else
            ! if the first dimension is the same, compare the second dimension
            if ((b_i - b_j) .lt. -1e-6) then
                smallerThan = .true.
            else if ((b_i - b_j) .gt. 1e-6) then
                smallerThan = .false.
            else
                ! if the second dimension is the same, compare the third dimension
                if ((c_i - c_j) .lt. -1e-6) then
                    smallerThan = .true.
                else
                    smallerThan = .false.
                end if
            end if
        end if

        return
    end function

    ! Calculate the distance squared between two vectors
    ! vector1, vector2: two vectors being calculated
    ! max_dist_square: if the distance exceeds this value, no need to further calculate
    ! dist_square: return the distance if it does not exceed the max_dist_square
    function vector_dist_square(vector1, vector2, max_dist_square) result(dist_square)
        real, dimension(3), intent(in) :: vector1, vector2
        real :: max_dist_square, dist_square, vector
        integer :: i

        dist_square = 0.0E0_dp
        do i = 1, 3
            vector = vector1(i) - vector2(i) 
            dist_square = dist_square + vector * vector
            if (dist_square .gt. max_dist_square) return
        end do
    end function vector_dist_square

    ! Calculate the distance squared between one point and the cube of a node
    ! mol_tree: the pointer the tree
    ! vector: the coordinate
    ! current_node: the pointer to the node whose cube is being used in the calculation
    ! max_dist_square: if the distance exceeds this value, no need to further calculate
    ! dist_square: return the distance if it does not exceed the max_dist_square
    function dist_square_from_point_to_cube(mol_tree, vector, current_node, max_dist_square) result(dist_square)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node 
        real, dimension(3), intent(in) :: vector
        real :: dist_square, max_dist_square, upper, lower
        integer :: i_dim, cut_dim

        dist_square = 0.0E0_dp
        do i_dim = 1, 3
            if (i_dim .eq. 1) then
                cut_dim = cut_dim_prev(current_node%cut_dim)
                lower = current_node%cube%lower
                upper = current_node%cube%upper
            else if (i_dim .eq. 2) then
                cut_dim = cut_dim_prev(cut_dim) 
                if (associated(current_node%parent_node)) then
                    lower = current_node%parent_node%cube%lower
                    upper = current_node%parent_node%cube%upper
                else
                    lower = mol_tree%cube(cut_dim)%lower
                    upper = mol_tree%cube(cut_dim)%upper
                end if
            else
                cut_dim = cut_dim_prev(cut_dim)
                if (current_node%height .ge. 3) then
                    lower = current_node%parent_node%parent_node%cube%lower
                    upper = current_node%parent_node%parent_node%cube%upper
                else
                    lower = mol_tree%cube(cut_dim)%lower
                    upper = mol_tree%cube(cut_dim)%upper
                end if
            end if
    
            if (vector(cut_dim) .lt. lower) then
                dist_square = dist_square + (vector(cut_dim) - lower)**2
            else if (vector(cut_dim) .gt. upper) then
                dist_square = dist_square + (vector(cut_dim) - upper)**2
            end if
            
            if (dist_square .gt. max_dist_square) return
        end do

    end function

    ! Calculate the next dimension
    ! cut_dim: the cutting dimension
    ! next_cut_dim: the next cutting dimension
    function cut_dim_next(cut_dim) result(next_cut_dim)
        integer, intent(in) :: cut_dim
        integer :: next_cut_dim
        
        if (cut_dim .eq. 1) then
            next_cut_dim = 2
        else if (cut_dim .eq. 2) then
            next_cut_dim = 3
        else
            next_cut_dim = 1
        end if
    end function

    ! Return the previous dimension
    ! cut_dim: the cutting dimension
    ! previous_cut_dim: the previous cutting dimension
    function cut_dim_prev(cut_dim) result(previous_cut_dim)
        integer, intent(in) :: cut_dim
        integer :: previous_cut_dim

        if (cut_dim .eq. 1) then
            previous_cut_dim = 3
        else if (cut_dim .eq. 2) then
            previous_cut_dim = 1
        else
            previous_cut_dim = 2
        end if
    end function cut_dim_prev
        
    ! Determine whether to go left or right from a parent node
    ! coord_to_add: 3d real array, the coordinate of the node to be added/searched/deleted
    ! coord_to_compare: 3d real array, the coodinate of the node to be compared (typically the current node)
    ! cut_dim: the cutting dimension of the current node
    ! lLeft: return true if it should go left, false if it should go right
    ! this fxn is similar to smaller_than(), but the input array structure is different
    ! for efficient reason, we use l_left in the node search, insertion, deletion and range_search
    function l_left(coord_to_add, coord_to_compare, cut_dim) result(lLeft)
        real, dimension(3) :: coord_to_add, coord_to_compare
        integer :: cut_dim, next_cut_dim, next_next_cut_dim
        logical :: lLeft

        if ((coord_to_add(cut_dim) - coord_to_compare(cut_dim)) .lt. -1e-6) then
            lLeft = .true.
        else if ((coord_to_add(cut_dim) - coord_to_compare(cut_dim)) .gt. 1e-6) then
            lLeft = .false.
        else
            next_cut_dim = cut_dim_next(cut_dim)
            if ((coord_to_add(next_cut_dim) - coord_to_compare(next_cut_dim)) .lt. -1e-6) then
                lLeft = .true.
            else if ((coord_to_add(next_cut_dim) - coord_to_compare(next_cut_dim)) .gt. 1e-6) then
                lLeft = .false.
            else
                next_next_cut_dim = cut_dim_next(next_cut_dim)
                if ((coord_to_add(next_next_cut_dim) - coord_to_compare(next_next_cut_dim)) .lt. -1e-6) then
                    lLeft = .true.
                else
                    lLeft = .false.
                end if
            end if
        end if
    end function l_left


    ! Set the cube for current_node
    ! mol_tree: the pointer to the tree
    ! current_node: the pointer to the current node
    ! lLeft: whether the current node is the left or right child of its parent
    subroutine set_cube(mol_tree, current_node, lLeft)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        integer :: cut_dim
        logical, intent(in) :: lLeft
        real :: lower, upper

        current_node%l_cube_updated = .true.

        ! if it's the root, the cube is for z-dimension
        if (.not. associated(current_node%parent_node)) then
            current_node%cube%lower = mol_tree%cube(3)%lower
            current_node%cube%upper = mol_tree%cube(3)%upper
            return
        end if

        cut_dim = current_node%parent_node%cut_dim !< cut_dim is the parent's cut_dim

        ! Determine whether there is great grandparent for the current node
        if (current_node%height .ge. 4) then
            lower = current_node%parent_node%parent_node%parent_node%cube%lower
            upper = current_node%parent_node%parent_node%parent_node%cube%upper
        else
            lower = mol_tree%cube(cut_dim)%lower
            upper = mol_tree%cube(cut_dim)%upper
        end if

        if (lLeft) then
            current_node%cube%lower = lower
            current_node%cube%upper = current_node%parent_node%coord(cut_dim)
        else    
            current_node%cube%lower = current_node%parent_node%coord(cut_dim)
            current_node%cube%upper = upper 
        end if
        
    end subroutine set_cube

    ! to initialize a kd-tree
    ! mol_tree: the pointer to the tree being calculated
    ! ibox: the box that this tree represents
    ! cubeP: the cube of the root
    ! mol_tree_initialized: return the pointer to the tree being initialized
    function init_tree(mol_tree, ibox, cubeP) result(mol_tree_initialized)
        type(tree), pointer :: mol_tree, mol_tree_initialized
        integer :: ibox
        type(interval), pointer :: cubeP(:)
        if (associated(mol_tree)) call err_exit(__FILE__,__LINE__,'kd-tree already existed!',1)
        allocate(mol_tree)
        allocate(mol_tree%cube(3))
        allocate(mol_tree%bound(3))
        allocate(mol_tree%bound_all(3))
        mol_tree%cube = cubeP
        mol_tree%node_num = 0
        mol_tree%height = 0
        mol_tree%box = ibox
        nullify(mol_tree%tree_root)
        mol_tree_initialized => mol_tree
    end function init_tree

    ! empty all the nodes in the tree
    subroutine empty_tree(mol_tree)
        type(tree), pointer :: mol_tree
        if (associated(mol_tree)) call empty_node(mol_tree%tree_root)
        deallocate(mol_tree%cube)
        deallocate(mol_tree%bound)
        deallocate(mol_tree%bound_all)
        deallocate(mol_tree)
    end subroutine empty_tree

    ! recursively empty the nodes
    recursive subroutine empty_node(current_node)
        type(tree_node), pointer :: current_node
        if (associated(current_node%left_node)) then
            call empty_node(current_node%left_node)
            nullify(current_node%left_node)
        end if

        if (associated(current_node%right_node)) then
            call empty_node(current_node%right_node)
            nullify(current_node%right_node)
        end if

        if (associated(current_node%cube)) deallocate(current_node%cube) 
        deallocate(current_node) 
    end subroutine empty_node
   
    ! Insert a node to the tree
    ! mol_tree: the pointer to the tree where the node is being inserted
    ! coord_to_add: 3d real array, the coordinate of the node to be added
    ! ichain: the chain number of the molecule
    ! ibead: the bead number of the molecule 
    ! ix, iy, iz: which periodic image it is on the x, y or z axis, a value between -1, 1 and 0
    ! mol_tree_inserted: return the pointer to the tree after the insertion
    function insert_node(mol_tree, coord_to_add, ichain, ibead, ix, iy, iz) result(mol_tree_inserted)
        type(tree), pointer :: mol_tree, mol_tree_inserted
        real, dimension(3) :: coord_to_add
        integer, intent(in) :: ichain, ibead, ix, iy, iz
        integer :: i
    
        if (.not. associated(mol_tree)) call err_exit(__FILE__,__LINE__,'kd-tree has not been initialized yet!',1)

        mol_tree%tree_root => insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead, ix, iy, iz, &
                                mol_tree%tree_root, 1, 1, .true., null())

        ! Update the boundary
        do i = 1, 3
            if (coord_to_add(i) .gt. mol_tree%bound_all(i)%upper) mol_tree%bound_all(i)%upper = coord_to_add(i) + 0.01
            if (coord_to_add(i) .lt. mol_tree%bound_all(i)%lower) mol_tree%bound_all(i)%lower = coord_to_add(i) - 0.01
        end do

        mol_tree_inserted => mol_tree

    end function insert_node

    ! Insert the coordinates coord_to_add to the proper node
    ! mol_tree: pointer to the tree
    ! coord_to_add: 3d real array, the coordinate of the node to be added
    ! ichain: the chain number of the molecule
    ! ibead: the bead number of the molecule
    ! ix, iy, iz: which periodic image it is on the x, y or z axis, a value between -1, 1 and 0 
    ! current_node: the current searching node
    ! height: the height of the current node
    ! cut_dim: the cutting dimension of the current node
    ! lLeft: true if it's the left of its parent node
    ! parent_node: the pointer to the parent node
    ! node_added: return the current searching node
    recursive function insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead, ix, iy, iz, current_node, height&
                        ,cut_dim, lLeft, parent_node) result(node_added)
        type(tree), pointer :: mol_tree
        real, dimension(3) :: coord_to_add
        integer :: ichain, ibead, ix, iy, iz
        type(tree_node), pointer :: current_node, parent_node, node_added
        integer :: height
        integer, intent(in) :: cut_dim
        logical, intent(in) :: lLeft
        integer :: next_cut_dim
        next_cut_dim = cut_dim_next(cut_dim)

        if (.not. associated(current_node)) then
            allocate(current_node)
            current_node%coord = coord_to_add
            current_node%ichain = ichain
            current_node%ibead = ibead
            current_node%ix = ix
            current_node%iy = iy
            current_node%iz = iz
            current_node%l_cube_updated = .false.
            current_node%parent_node => parent_node
            current_node%cut_dim = cut_dim
            current_node%height = height
            allocate(current_node%cube)
            if (.not. associated(current_node%parent_node)) call set_cube(mol_tree, current_node, lLeft)
            mol_tree%node_num = mol_tree%node_num + 1
            if (height .gt. mol_tree%height) mol_tree%height = height
            nullify(current_node%left_node, current_node%right_node)
        else if ((vector_dist_square(current_node%coord, coord_to_add, 1e-6) .lt. 1e-6) &
                .and. (current_node%ichain .eq. ichain) .and. (current_node%ibead .eq. ibead) &
                .and. (current_node%ix .eq. ix) .and. (current_node%iy .eq. iy) .and. (current_node%iz .eq. iz)) then
            write(*,*) "current_node: coordinates", current_node%coord
            write(*,*) "current_node: ichain, ibead", current_node%ichain, current_node%ibead
            write(*,*) "current_node: ix, iy, iz", current_node%ix, current_node%iy, current_node%iz
            write(*,*) "coord_to_add: coordinates", coord_to_add
            write(*,*) "coord_to_add: ichain, ibead", ichain, ibead
            write(*,*) "coord_to_add: ix, iy, iz", ix, iy, iz
            call err_exit(__FILE__,__LINE__,'Duplicate keys added when constructing the kd-tree',1)
        else if (l_left(coord_to_add, current_node%coord, cut_dim)) then 
            current_node%left_node => insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead &
                , ix, iy, iz, current_node%left_node, height+1, next_cut_dim, .true., current_node)
        else
            current_node%right_node => insert_node_into_tree(mol_tree, coord_to_add, ichain, ibead &
                , ix, iy, iz, current_node%right_node, height+1, next_cut_dim, .false., current_node)
        end if

        node_added => current_node

    end function insert_node_into_tree
    
    ! Return the minimum node in the tree at dimension min_dim
    ! current_node: the current node being examined
    ! min_dim: the dimension where the minimum is
    ! cut_dim: the current dimension
    ! min_node: return the pointer to the minimum node
    recursive function find_min(current_node, min_dim, cut_dim) result(min_node)
        type(tree_node), pointer :: current_node, min_node, min_nodeLeft, min_node_right
        integer, intent(in) :: min_dim, cut_dim
        real :: min_node_coord(3)
        logical :: lLeft, lRight
        
        if (.not. associated(current_node)) min_node => null()
        
        if (cut_dim .eq. min_dim) then
            if (.not. associated(current_node%left_node)) then 
                min_node => current_node
            else 
                min_node => find_min(current_node%left_node, min_dim, cut_dim_next(cut_dim))
            end if
        else
            min_node => current_node
            min_node_coord = current_node%coord
            lLeft = associated(current_node%left_node)
            lRight = associated(current_node%right_node)

            ! if there is sub-tree, compare the value with left and right sub-tree minimums
            if (lLeft) then
                min_nodeLeft => find_min(current_node%left_node, min_dim, cut_dim_next(cut_dim))
                if (l_left(min_nodeLeft%coord, min_node_coord, min_dim)) then
                    min_node => min_nodeLeft
                    min_node_coord = min_node%coord
                end if
            end if

            if (lRight) then
                min_node_right => find_min(current_node%right_node, min_dim, cut_dim_next(cut_dim))
                if (l_left(min_node_right%coord, min_node_coord, min_dim)) min_node => min_node_right
            end if
        end if

    end function find_min

    ! Return the maximum node in the tree at dimension min_dim
    ! current_node: the current node being examined
    ! max_dim: the dimension where the maximum is
    ! cut_dim: the current dimension
    ! max_node: return the pointer to the maximum node
    recursive function find_max(current_node, max_dim, cut_dim) result(max_node)
        type(tree_node), pointer :: current_node, max_node, max_node_left, max_nodeRight
        integer, intent(in) :: max_dim, cut_dim
        real :: max_nodeCoord(3)
        logical :: lLeft, lRight
        
        if (.not. associated(current_node)) max_node => null()
        
        if (cut_dim .eq. max_dim) then
            if (.not. associated(current_node%right_node)) then 
                max_node => current_node
            else 
                max_node => find_max(current_node%right_node, max_dim, cut_dim_next(cut_dim))
            end if
        else
            max_node => current_node
            max_nodeCoord = current_node%coord
            lLeft = associated(current_node%left_node)
            lRight = associated(current_node%right_node)

            ! if there is sub-tree, compare the value with left and right sub-tree minimums
            if (lLeft) then
                max_node_left => find_max(current_node%left_node, max_dim, cut_dim_next(cut_dim))
                if (.not. l_left(max_node_left%coord, max_nodeCoord, max_dim)) then
                    max_node => max_node_left
                    max_nodeCoord = max_node%coord
                end if
            end if

            if (lRight) then
                max_nodeRight => find_max(current_node%right_node, max_dim, cut_dim_next(cut_dim))
                if (.not. l_left(max_nodeRight%coord, max_nodeCoord, max_dim)) max_node => max_nodeRight
            end if
        end if

    end function find_max

    ! Update the height of the tree
    ! mol_tree: the pointer to the tree
    subroutine update_tree_height(mol_tree)
        type(tree), pointer :: mol_tree
    
        mol_tree%height = 0
        if (.not. associated(mol_tree)) return
        call update_height_in_tree(mol_tree, mol_tree%tree_root, 1, 1)
    end subroutine update_tree_height

    ! Recursively find the height of the tree
    ! mol_tree: the pointer to the tree
    ! current_node: the current node being worked on
    ! cut_dim: the cutting dimension of the current node
    ! height: the height of the current node
    recursive subroutine update_height_in_tree(mol_tree, current_node, cut_dim, height)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        integer :: cut_dim, next_cut_dim, height

        next_cut_dim = cut_dim_next(cut_dim)

        if (height .gt. mol_tree%height) mol_tree%height = height
        if (associated(current_node%left_node)) call update_height_in_tree(mol_tree, current_node%left_node, next_cut_dim, height+1)
        if (associated(current_node%right_node)) call update_height_in_tree(mol_tree, current_node%right_node, next_cut_dim, height+1)
         
    end subroutine update_height_in_tree

    ! Delete a node in the tree, return the pointer to the tree
    ! mol_tree: the pointer to the tree
    ! coord_to_delete: 3d coordinates of the nodes to be deleted
    ! mol_tree_deleted: return the pointer to the tree
    function delete_node(mol_tree, coord_to_delete) result(mol_tree_deleted)
        type(tree), pointer :: mol_tree, mol_tree_deleted
        real, dimension(3), intent(in) :: coord_to_delete
        real, dimension(3) :: coord_to_compare
        type(tree_node), pointer :: node_replaced

        coord_replaced = [10000, 10000, 10000] !< set to unrealistic values       
        coord_to_compare = coord_replaced

        mol_tree%tree_root => delete_node_in_tree(mol_tree, mol_tree%tree_root, coord_to_delete, 1, .false.)      
 
        mol_tree%node_num = mol_tree%node_num - 1
        mol_tree_deleted => mol_tree
      
        if (vector_dist_square(coord_to_compare, coord_replaced, 1e-6) .gt. 1e-6) then
            node_replaced => search_node_in_tree(mol_tree, mol_tree%tree_root, coord_replaced, 1)
            node_replaced => reset_cube_status(node_replaced)
        end if

    end function delete_node

    ! Recursively delete a node in the tree
    ! mol_tree: the pointer to the tree
    ! current_node: the node currently being examined
    ! coord_to_delete: 3d coordinates of the nodes to be deleted
    ! cut_dim: the cutting dimension of the current node
    ! l_update_cube: whether to update cube, false until the the first occurence of node deletion
    ! node_deleted: return the pointer to the node AFTER the deletion
    recursive function delete_node_in_tree(mol_tree, current_node, coord_to_delete, cut_dim, l_update_cube) result(node_deleted)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node, node_deleted, node_replaced
        real, dimension(3), intent(in) :: coord_to_delete
        integer, intent(in) :: cut_dim
        logical, intent(in) :: l_update_cube
        integer :: next_cut_dim

        next_cut_dim = cut_dim_next(cut_dim)

        if (.not. associated(current_node)) then
            write(*,*) coord_to_delete
            call err_exit(__FILE__,__LINE__,'Did not find the node!',1)
        end if

        if ((abs(current_node%coord(1) - coord_to_delete(1)) .lt. 1e-6) .and. &
            (abs(current_node%coord(2) - coord_to_delete(2)) .lt. 1e-6) .and. &
            (abs(current_node%coord(3) - coord_to_delete(3)) .lt. 1e-6)) then

            ! if this is the node to delete
            if (associated(current_node%right_node)) then
                ! if there's only right node, use the minimum from the right sub-tree
                node_replaced => find_min(current_node%right_node, cut_dim, next_cut_dim)

                current_node%coord = node_replaced%coord
                current_node%ichain = node_replaced%ichain
                current_node%ibead = node_replaced%ibead
                current_node%ix = node_replaced%ix
                current_node%iy = node_replaced%iy
                current_node%iz = node_replaced%iz
                
                ! record the coordinates to coord_replaced if l_update_cube == .false.
                if (.not. l_update_cube) coord_replaced = node_replaced%coord

                current_node%right_node => delete_node_in_tree(mol_tree, current_node%right_node, current_node%coord, &
                                        next_cut_dim, .true.)
            else if (associated(current_node%left_node)) then
                ! if there's only left node, swap sub-trees, and use the minimum from the next right-tree
                node_replaced => find_min(current_node%left_node, cut_dim, next_cut_dim)
                current_node%coord = node_replaced%coord
                current_node%ichain = node_replaced%ichain
                current_node%ibead = node_replaced%ibead
                current_node%ix = node_replaced%ix
                current_node%iy = node_replaced%iy
                current_node%iz = node_replaced%iz
                
                ! record the coordinates to coord_replaced if l_update_cube == .false.
                if (.not. l_update_cube) coord_replaced = node_replaced%coord

                current_node%right_node => delete_node_in_tree(mol_tree, current_node%left_node, current_node%coord, &
                                        next_cut_dim, .true.)
                current_node%left_node => null()
            else
                ! if it's only a leaf, just remove
                deallocate(current_node%cube)
                deallocate(current_node)
                node_deleted => null()
                return
            end if

        ! if this is not the node to delete, search for it
        else if (l_left(coord_to_delete,  current_node%coord, cut_dim)) then
            current_node%left_node => delete_node_in_tree(mol_tree, current_node%left_node, coord_to_delete, next_cut_dim, l_update_cube)
        else
            current_node%right_node => delete_node_in_tree(mol_tree, current_node%right_node, coord_to_delete, next_cut_dim, l_update_cube)
        end if

        ! Return the node 
        node_deleted => current_node

    end function delete_node_in_tree

    ! Reset all the l_cube_updated to be .false. recursively
    ! current_node: the starting node to be updated
    ! node_updated: return the updated node pointer
    recursive function reset_cube_status(current_node) result(node_updated)
        type(tree_node), pointer :: current_node, node_updated

        current_node%l_cube_updated = .false.
        if (associated(current_node%left_node)) current_node%left_node => reset_cube_status(current_node%left_node)        
        if (associated(current_node%right_node)) current_node%right_node => reset_cube_status(current_node%right_node)        

        node_updated => current_node
    end function reset_cube_status

    ! Update the cube for all the sub-trees of coord_to_update
    ! mol_tree: the pointer to the tree
    ! coord_to_update: 3d coordinates of the node whose sub-trees need to be updated
    ! mol_treeUpdated: return the pointer to the tree after the cube update
    function update_cube(mol_tree, coord_to_update) result(mol_treeUpdated)
        type(tree), pointer :: mol_tree, mol_treeUpdated
        real, dimension(3) :: coord_to_update
        type(tree_node), pointer :: node_to_update
        
        node_to_update => search_node_in_tree(mol_tree, mol_tree%tree_root, coord_to_update, 1)
        node_to_update => update_cube_in_tree(mol_tree, node_to_update)
        mol_treeUpdated => mol_tree
    end function update_cube

    ! Update the cube for all the sub-trees of coord_to_update
    ! mol_tree: the pointer to the tree
    ! current_node: the node currently working on
    ! node_updated: the updated current node
    recursive function update_cube_in_tree(mol_tree, current_node) result(node_updated)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node, node_updated

        if (associated(current_node%left_node)) then
            call set_cube(mol_tree, current_node%left_node, .true.)
            current_node%left_node => update_cube_in_tree(mol_tree, current_node%left_node)
        end if

        if (associated(current_node%right_node)) then
            call set_cube(mol_tree, current_node%right_node, .false.)
            current_node%right_node => update_cube_in_tree(mol_tree, current_node%right_node)
        end if

        node_updated => current_node

    end function update_cube_in_tree

    ! Search a node in the tree
    ! current_node: the node currently working on
    ! coord_to_search: 3d coordinates of the node to be searched
    ! cut_dim: the cutting dimension of the current_node
    ! node_found: the node found
    recursive function search_node_in_tree(mol_tree, current_node, coord_to_search, cut_dim) result(node_found)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node, node_found
        real, dimension(3) :: coord_to_search
        integer, intent(in) :: cut_dim

        if (.not. associated(current_node)) then
            node_found => null()
            return
        end if

        if ((abs(current_node%coord(1) - coord_to_search(1)) .lt. 1e-6) .and. &
            (abs(current_node%coord(2) - coord_to_search(2)) .lt. 1e-6) .and. &
            (abs(current_node%coord(3) - coord_to_search(3)) .lt. 1e-6)) then

            node_found => current_node

        else if (l_left(coord_to_search, current_node%coord, cut_dim)) then
            node_found => search_node_in_tree(mol_tree, current_node%left_node, coord_to_search, cut_dim_next(cut_dim))
        else    
            node_found => search_node_in_tree(mol_tree, current_node%right_node, coord_to_search, cut_dim_next(cut_dim))
        end if
        
    end function search_node_in_tree

    ! Check the validity of the tree, for the coordinates only
    ! mol_tree: the pointer to the tree
    ! lcutcm: whether it is COM-kdtree
    ! lValid: whether the tree structure is correct
    function check_tree_coord(mol_tree, lcutcm) result(lValid)
        type(tree), pointer :: mol_tree
        logical :: lValid, lcutcm
        integer :: ibox, i

        ibox = mol_tree%box

        lValid = check_coord_in_tree(mol_tree%tree_root, 1, lcutcm, ibox)
        if (.not. lValid) then
            do i = 1, 3
                write(*,*) "mol_tree%bound_all", mol_tree%bound_all(i)%lower, mol_tree%bound_all(i)%upper
            end do
        end if
    end function check_tree_coord
    
    ! Check the validity of the node in the tree
    ! current_node: the node being examined
    ! lcutcm: whether it is COM-kdtree
    ! ibox: the box of the tree
    ! lValid: whether the tree structure at this node is correct
    recursive function check_coord_in_tree(current_node, cut_dim, lcutcm, ibox) result(lValid)
        use sim_system, only : xcm, ycm, zcm, boxlx, boxly, boxlz
        type(tree_node), pointer :: current_node, max_node_left, min_node_right
        logical :: lValid, lcutcm
        integer :: cut_dim, next_cut_dim, ibox, ichain
        real :: x_coord, y_coord, z_coord, xcm_coord, ycm_coord, zcm_coord

        ! if the cutting dimension does not match
        if (current_node%cut_dim .ne. cut_dim) then
            lValid = .false.
            return
        end if

        ! check whether the coordinate matches with the coordinates in the unordered array
        if (lcutcm) then
            ichain = current_node%ichain
            x_coord = current_node%coord(1)
            y_coord = current_node%coord(2)
            z_coord = current_node%coord(3)
            xcm_coord = xcm(ichain) + boxlx(ibox) * current_node%ix
            ycm_coord = ycm(ichain) + boxly(ibox) * current_node%iy
            zcm_coord = zcm(ichain) + boxlz(ibox) * current_node%iz

            if ((abs(x_coord-xcm_coord) .gt. 1e-6 ) .or. &
                (abs(y_coord-ycm_coord) .gt. 1e-6 ) .or. &
                (abs(z_coord-zcm_coord) .gt. 1e-6 )) then

                lValid = .false.
                return
            end if
        end if

        next_cut_dim = cut_dim_next(cut_dim)
        lValid = .true.

        if (associated(current_node%left_node)) then
            max_node_left => find_max(current_node%left_node, cut_dim, next_cut_dim)
            if (l_left(current_node%coord, max_node_left%coord, cut_dim)) then
                lValid = .false.
                return
            else
                lValid = check_coord_in_tree(current_node%left_node, next_cut_dim, lcutcm, ibox)
                if (.not. lValid) return
            end if
        end if
    
        if (associated(current_node%right_node)) then
            min_node_right => find_min(current_node%right_node, cut_dim, next_cut_dim)
            if ((l_left(min_node_right%coord, current_node%coord, cut_dim)) &
              .and. (.not. l_left(current_node%coord, min_node_right%coord, cut_dim))) then
                lValid = .false.
                return
            else
                lValid = check_coord_in_tree(current_node%right_node, next_cut_dim, lcutcm, ibox)
                if (.not. lValid) return
            end if
        end if
    end function check_coord_in_tree

    ! Check the validity of the cube in the tree
    ! mol_tree: the pointer to the tree
    ! lValid: whetehr the tree cube is correct
    function check_tree_cube(mol_tree) result(lValid)
        type(tree), pointer :: mol_tree
        logical :: lValid
        lValid = check_cube_in_tree(mol_tree, mol_tree%tree_root)
    end function check_tree_cube

    ! Check the validity of the cube in the current_node
    ! mol_tree: the pointer to the tree
    ! current_node: the node being examined
    ! lValid: whether the tree cube at this node is correct
    recursive function check_cube_in_tree(mol_tree, current_node) result(lValid)
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        logical :: lValid
        integer :: previous_cut_dim
        real :: upper, lower

        ! if it's the root
        if (.not. associated(current_node%parent_node)) then
            if ((current_node%cube%upper .eq. mol_tree%cube(3)%upper) &
                    .and. (current_node%cube%lower .eq. mol_tree%cube(3)%lower)) then
                lValid = .true.
            else
                lValid = .false.
                return
            end if
 
            !if (.not. lValid) return
            if (associated(current_node%left_node)) lValid = check_cube_in_tree(mol_tree, current_node%left_node)
            if (.not. lValid) return
            if (associated(current_node%right_node)) lValid = check_cube_in_tree(mol_tree, current_node%right_node)
            return
        end if

        previous_cut_dim = current_node%parent_node%cut_dim

        ! if it's NOT the root, check the cube
        lower = current_node%cube%lower
        upper = current_node%cube%upper

        if (l_left(current_node%coord, current_node%parent_node%coord, previous_cut_dim)) then
            call set_cube(mol_tree, current_node, .true.)
        else
            call set_cube(mol_tree, current_node, .false.)
        end if

        if ((lower .ne. current_node%cube%lower) .or. (upper .ne. current_node%cube%upper)) then
            lValid = .false.
            return
        else
            lValid = .true.
        end if       

        ! recursively check its leaves
        if (associated(current_node%left_node)) lValid = check_cube_in_tree(mol_tree, current_node%left_node)
        if (.not. lValid) return
        if (associated(current_node%right_node)) lValid = check_cube_in_tree(mol_tree, current_node%right_node)
        return
        
    end function check_cube_in_tree

    ! Search for the coordinates within rcut with respect to a reference point ref_coord
    ! Return loverlap = true if there is an r larger than rmin
    ! otherwise return an array of eligible search_output
    ! mol_tree: pointer the tree
    ! ref_coord: the reference coordinate, whose interactions with all other beads are being calculated
    ! max_dim: maximum dimension of the search_output
    ! rmin: minimum distance allowed
    ! rcut: cutoff distance
    ! ichain: ref_coord is for ichain; do not search for bead in ichain
    ! ibead: ref_coord is for ibead in ichain; used to determine what type of bead it is and whether it has LJ or QQ interaction
    ! lsumup: true if called from sumup suborutine, accout for only interactions when chain number is greater than ichain
    !         false if called from energy or boltz subroutine, account for all the interactions when bead is not in ichain
    ! loverlap: return true if there's an overlap
    ! actual_dim: actual dimension of the search_output
    ! search_output_array: 3 or 6 * actual_dim real array which have information about the qualified beads to calculate interactions
    ! dist_calc_num: number of distance calculation used
    ! lPressure: whether called from pressure calculation, if so, have additional dimension of output arrays for r*uij
    subroutine range_search(mol_tree, ref_coord, max_dim, rmin, rcut, ichain, ibead, lsumup, &
        loverlap, actual_dim, search_output, dist_calc_num, lPressure)
        use sim_system, only : lcutcm
        type(tree), pointer :: mol_tree
        real, dimension(3), intent(in) :: ref_coord
        integer, intent(in) :: max_dim, ichain, ibead 
        real, intent(in) :: rmin, rcut
        logical, intent(in) :: lPressure
        logical :: lsumup, loverlap
        real, allocatable :: search_output(:, :)
        integer :: i, dist_calc_num, search_output_dim, actual_dim
        integer ::  i_search_output

        ! Initialize the variables
        if (allocated(search_output)) deallocate(search_output)

        ! If lcutcm == .false. (all beads are stored)
        ! dist_square, ichain, ibead for 3 dimensions, if lPressure, additional 3d for rxuij, ryuij and rzuij
        ! else if lcutcm == .true. (COMs are stored)
        ! 1 dimension (search_output ~= lij2), logical,
        ! indicating whether we should include the interaction between i and j
        if (lcutcm) then
            search_output_dim = 1
        else
            if (lPressure) then
                search_output_dim = 6
            else
                search_output_dim = 3
            end if
        end if

        allocate(search_output(search_output_dim, max_dim))

        if (lcutcm) search_output(1,:) = 0.0 !< 0.0 = false, 1.0 = true

        loverlap = .false.
        rmin_square = rmin * rmin
        rcut_square = rcut * rcut
        dist_calc_num = 0

        ! Start searching
        i_search_output = 1
        call range_search_in_tree(mol_tree, mol_tree%tree_root, ref_coord, 1, ichain, ibead &
            , lsumup, loverlap, dist_calc_num, lPressure, search_output, i_search_output)

        actual_dim = i_search_output - 1

    end subroutine range_search

    ! recursively do the range search
    ! mol_tree: pointer the tree
    ! current_node: the node currently working on
    ! ref_coord: the reference coordinate, whose interactions with all other beads are being calculated
    ! cut_dim: the cutting dimension of the current node
    ! ichain: ref_coord is for ichain; do not search for bead in ichain
    ! ibead: ref_coord is for ibead in ichain; used to determine what type of bead it is and whether it has LJ or QQ interaction
    ! lsumup: true if called from sumup suborutine, accout for only interactions when chain number is greater than ichain
    !         false if called from energy or boltz subroutine, account for all the interactions when bead is not in ichain
    ! loverlap: return true if there's an overlap
    ! dist_calc_num: number of distance calculation used
    ! lPressure: whether called from pressure calculations, add additional 3d output for r*uij
    ! search_output: 3 or 6 * max_dim array which has information about the qualified beads to calculate interactions
    ! i_search_output: the current index of search_output
    recursive subroutine range_search_in_tree(mol_tree, current_node, ref_coord, cut_dim, ichain, ibead, lsumup&
                , loverlap, dist_calc_num, lPressure, search_output, i_search_output)
        use sim_system,only:io_output,lcutcm

        type(tree), pointer :: mol_tree
        real, dimension(3), intent(in) :: ref_coord
        integer, intent(in) :: cut_dim, ichain, ibead
        type(tree_node), pointer :: current_node, node_closer, node_farther
        logical, intent(in) :: lPressure
        real :: search_output(:, :), dist_square
        integer :: i_search_output, dist_calc_num, next_cut_dim
        logical :: lsumup, loverlap, lLeft

        ! if overlap, then return
        if (loverlap) return

        ! calculate the distance and include this point if it's within rcut
        ! only for the intermolecular part
        ! if .not. lsumup, account for chain number different from ichain
        ! if lsumup, only account for chain number greater than ichain
        if ((.not. lsumup .and. current_node%ichain .ne. ichain) .or. (lsumup .and. current_node%ichain .gt. ichain)) then
            dist_square = vector_dist_square(current_node%coord, ref_coord, rcut_square)
            dist_calc_num = dist_calc_num + 1
            if (dist_square .lt. rmin_square) then
                ! It is debatable whether we need to call overlap if two beads that are closer than rmin
                ! but they don't have interactions (e.g. H and O in TIP4P water)
                !if (check_interaction(ichain, ibead, current_node%ichain, current_node%ibead)) then
                loverlap = .true.
                !    return
                !end if
            else if (dist_square .le. rcut_square) then
                if (lcutcm) then
                    search_output(1, current_node%ichain) = 1.0 !< 1.0 = true, 0.0 = false
                else
                    search_output(1, i_search_output) = dist_square
                    search_output(2, i_search_output) = current_node%ichain
                    search_output(3, i_search_output) = current_node%ibead

                    if (lPressure) then
                        search_output(4, i_search_output) = ref_coord(1) - current_node%coord(1)
                        search_output(5, i_search_output) = ref_coord(2) - current_node%coord(2)
                        search_output(6, i_search_output) = ref_coord(3) - current_node%coord(3)
                    end if

                    i_search_output = i_search_output + 1
                end if

            end if
        end if

        ! recursively search left or right in a more promising order
        next_cut_dim = cut_dim_next(cut_dim)
       
        if (l_left(ref_coord, current_node%coord, cut_dim)) then 
            node_closer => current_node%left_node
            node_farther => current_node%right_node
            lLeft = .true.
        else
            node_closer => current_node%right_node
            node_farther => current_node%left_node
            lLeft = .false.
        end if

        if (associated(node_closer)) then
            if (.not. node_closer%l_cube_updated) call set_cube(mol_tree, node_closer, lLeft)
            call range_search_in_tree(mol_tree, node_closer, ref_coord, next_cut_dim, ichain, ibead, lsumup, loverlap&
                    , dist_calc_num, lPressure, search_output, i_search_output)
        end if

        ! search on the farther node if necessary
        if (associated(node_farther)) then
            if (.not. node_farther%l_cube_updated) call set_cube(mol_tree, node_farther, .not. lLeft)
            dist_square = dist_square_from_point_to_cube(mol_tree, ref_coord, node_farther, rcut_square)
            if (dist_square .le. rcut_square) then
                call range_search_in_tree(mol_tree, node_farther, ref_coord, next_cut_dim, ichain, ibead, lsumup, loverlap&
                        , dist_calc_num, lPressure, search_output, i_search_output)
            end if 
        end if
    end subroutine range_search_in_tree

    ! Find the k-th element of an unordered array; used to find the median of the array
    ! k: the k-th element to be found
    ! coord_array1, 2 and 3: array of coordinates, based on which the sorting is done;
    !                        first compare coordinates in array 1, if ties, compare those in array 2
    ! index_array: the array with all the indices of the array
    ! bead_index: the index of the k-th element found (its index in the index_array)
    subroutine find_k(k, coord_array1, coord_array2, coord_array3, index_array, bead_index)
        integer, intent(in) :: k
        real, intent(in out) :: coord_array1(:), coord_array2(:), coord_array3(:)
        integer, intent(in out) :: index_array(:)
        integer, intent(out) :: bead_index

        integer :: hi, lo, i, j

        lo = 1
        hi = size(index_array)

        do while (hi .gt. lo)

            call partition(coord_array1, coord_array2, coord_array3, index_array, lo, hi, j)

            if (j .gt. k) then
                ! all the elements whose indices are greater than j are greater than the k-th element
                ! recursively search in the left side
                hi = j - 1
            else if (j .lt. k) then
                ! all the elements whose indices are smaller than i are greater than the k-th element
                ! recursively search in the right side
                lo = j + 1
            else
                ! j .eq. k, which means we find the k-th element
                ! return the k-th (j-th) element
                bead_index = index_array(j)
                return
 
            end if
        end do
 
        bead_index = index_array(k)

        return
    end subroutine find_k

    ! subroutine used in the selection algorithm to find the k-th element
    ! coord_array1, 2 and 3: same as find_k
    ! index_array: same as find_k
    ! lo, hi: the low and high bound of the selection algorithm search
    ! j: the j-th element found
    subroutine partition(coord_array1, coord_array2, coord_array3, index_array, lo, hi, j)
        real, intent(in out) :: coord_array1(:), coord_array2(:), coord_array3(:)
        integer, intent(in out) :: index_array(:)
        integer, intent(in) :: lo, hi
        integer, intent(out) :: j

        integer :: i

        i = lo + 1
        j = hi

        do while (.true.)

            ! i loops from left to right until either reaches the end or find a value that is greater than a[lo]
            do while (smallerThan(coord_array1(i),coord_array1(lo),coord_array2(i),coord_array2(lo),coord_array3(i),coord_array3(lo)))
                i = i + 1
                if (i .ge. hi) exit
            end do

            ! j loops from right to left until either reaches the end or find a value that is smaller than a[lo]
            do while (smallerThan(coord_array1(lo),coord_array1(j),coord_array2(lo),coord_array2(j),coord_array3(lo),coord_array3(j)))
                j = j - 1
                if (j .le. lo) exit
            end do

            ! if i and j cross, exit the while loop
            if (i .ge. j) exit

            ! else, exchange the element of i and j and continue
            call swap_real_array(coord_array1, i, j)
            call swap_real_array(coord_array2, i, j)
            call swap_real_array(coord_array3, i, j)
            call swap_int_array(index_array, i, j)
        end do

        ! exchange the element of lo and j
        call swap_real_array(coord_array1, lo, j)
        call swap_real_array(coord_array2, lo, j)
        call swap_real_array(coord_array3, lo, j)
        call swap_int_array(index_array, lo, j)

        return
    end subroutine partition

    ! swap i and j element of a real array
    ! array: input array
    ! i, j: two indices where two elements need to be swapped
    subroutine swap_real_array(array, i, j)
        real :: array(:)
        integer :: i, j

        real :: real_temp

        real_temp = array(i)
        array(i) = array(j)
        array(j) = real_temp

        return
    end subroutine

    ! swap i and j element of an integer array
    ! same as above, except that this is for integer array
    subroutine swap_int_array(array, i, j)
        integer :: array(:)
        integer :: i, j

        integer :: int_temp

        int_temp = array(i)
        array(i) = array(j)
        array(j) = int_temp

        return
    end subroutine

    ! Initialize the tree coordinates by inserting the median coordinates of the cutting dimension
    ! mol_tree: the pointer to the tree
    ! rxu_tot, ryu_tot, rzu_tot: Arrays of x, y or z coordinates to be input
    ! index_array: N-dimensional array of indices, used to indicate/distinguish each bead
    ! tree_construction_list: the array of indices which indicate the sequence of inserting beads into the tree
    ! iList: the current index of tree_construction_list
    ! N_tot: the total dimension of the particles (originally), remains the same during recursive calls
    ! N: the number of total beads, also the "true" dimension of all the arrays
    ! cut_dim: the cutting dimension of the current node
    ! l_mpi: whether MPI is used to construct the coordinates
    ! i_numprocs: if MPI is used, the number of processors that has already been filled in with part of the list
    ! n_index_array: if MPI is used, the index_array for specific processor
    ! index_array_numprocs: if MPI is used, 2d-array (1st dim: index array index; 2nd dim: the number of cores) which
    !                       has information about the index array each processor works on
    ! i_level: if MPI is used, the cutting dimension of the bead that each processor starts with
    recursive subroutine init_tree_coordinates(rxu_tot, ryu_tot, rzu_tot, index_array, tree_construct_list, iList, N_tot, N, cut_dim &
            , l_mpi, i_numprocs, n_index_array, index_array_numprocs, i_level)
        use sim_system, only : numprocs, myid
        integer, intent(in) :: N_tot, N
        integer :: iList, bead_index
        real, intent(in) :: rxu_tot(:), ryu_tot(:), rzu_tot(:)
        integer, dimension(N_tot) :: tree_construct_list
        real, dimension(N) :: coord_sort1, coord_sort2, coord_sort3
        integer, dimension(N) :: index_array
        integer, allocatable :: index_array_left(:), index_array_right(:)
        integer :: cut_dim, next_cut_dim, median_index
        integer :: i, nLeft, nRight
        logical :: l_mpi !< if true, do not recursively process the array, stop after the array is sorted
        integer :: i_numprocs, i_level, num_array
        integer :: n_index_array(:), index_array_numprocs(:, :)

        if (cut_dim .eq. 1) then
            do i = 1, N
                bead_index = index_array(i)
                coord_sort1(i) = rxu_tot(bead_index)
                coord_sort2(i) = ryu_tot(bead_index)
                coord_sort3(i) = rzu_tot(bead_index)
            end do
            next_cut_dim = 2
        else if (cut_dim .eq. 2) then
            do i = 1, N
                bead_index = index_array(i)
                coord_sort1(i) = ryu_tot(bead_index)
                coord_sort2(i) = rzu_tot(bead_index)
                coord_sort3(i) = rxu_tot(bead_index)
            end do
            next_cut_dim = 3
        else
            do i = 1, N
                bead_index = index_array(i)
                coord_sort1(i) = rzu_tot(bead_index)
                coord_sort2(i) = rxu_tot(bead_index)
                coord_sort3(i) = ryu_tot(bead_index)
            end do
            next_cut_dim = 1
        end if 

        ! find and record the median node into tree_construct_list
        ! find median
        if (mod(N, 2) == 0) then
            median_index = N / 2
        else
            median_index = (N + 1) / 2
        end if

        call find_k(median_index, coord_sort1, coord_sort2, coord_sort3, index_array, bead_index)

        ! record the median in the tree_construct_list
        tree_construct_list(iList) = bead_index
        iList = iList + 1

        ! Count the number of elements in the left and right array
        nLeft = median_index - 1
        nRight = N - median_index

        ! if l_mpi is TRUE, it means that we are trying to divide the array and then distribute the sub-arrays
        ! to different cores for further recursive sorting
        ! of course only returns the sub-array when 2^(i_level) reaches numprocs
        if (l_mpi) then
            num_array = 2**(i_level)
            i_level = i_level + 1

            if (num_array .eq. numprocs) then
                n_index_array(i_numprocs) = nLeft
                index_array_numprocs(1: nLeft, i_numprocs) = index_array(1: nLeft)
                i_numprocs = i_numprocs + 1
                n_index_array(i_numprocs) = nRight
                index_array_numprocs(1: nRight, i_numprocs) = index_array(median_index + 1 : N)
                i_numprocs = i_numprocs + 1
                i_level = i_level - 1
                return
            end if
        end if

        ! Split arrays for the recursive call
        if (nLeft .ge. 1) then
            allocate(index_array_left(nLeft))
            index_array_left(1:nLeft) = index_array(1:nLeft)

            call init_tree_coordinates(rxu_tot, ryu_tot, rzu_tot, index_array_left, &
                tree_construct_list, iList, N_tot, nLeft, next_cut_dim, l_mpi, i_numprocs, n_index_array, index_array_numprocs, i_level)
        end if

        if (nRight .ge. 1) then
            allocate(index_array_right(nRight))
            index_array_right(1 : nRight) = index_array(median_index + 1 : N)

            call init_tree_coordinates(rxu_tot, ryu_tot, rzu_tot, index_array_right, &
                tree_construct_list, iList, N_tot, nRight, next_cut_dim, l_mpi, i_numprocs, n_index_array, index_array_numprocs, i_level)
        end if

        if (l_mpi) then
            i_level = i_level - 1
        end if

        return
         
    end subroutine init_tree_coordinates

    ! To construct/re-construct the kdtree
    !< ibox: for which simulation box
    !< iTree: if construct from scratch, iTree=ibox; if re-construct for volume move, iTree=nbox+1 or nbox+2
    !< lOutput: output the tree construction output if .true.
    subroutine construct_kdtree(ibox, iTree, lOutput)
        use sim_system
        use util_runtime,only:err_exit
        use util_mp,only:mp_set_displs,mp_allgather

        ! kd-tree variables
        logical :: lOutput, lAdd
        integer, intent(in) :: ibox, iTree
        type(interval), pointer :: cubeP(:)
        type(tree), pointer :: kd_tree
        real, allocatable :: rxu_tot(:), ryu_tot(:), rzu_tot(:)
        integer, allocatable :: bead_array(:), chain_array(:), ix_array(:), iy_array(:), iz_array(:), index_array(:), tree_construct_list(:)
        integer :: i, ix, iy, iz, ichain, ibead, N_max, N_tot, bead_index
        integer, allocatable :: nx(:), ny(:), nz(:)
        integer :: imolty, iList, nbead
        real :: xmin, xmax, ymin, ymax, zmin, zmax, xcoord, ycoord, zcoord, rbcut_plus_buffer
        real, dimension(3) :: coord

        ! MPI
        integer :: rcounts(numprocs),displs(numprocs),blocksize
        integer :: my_iList, my_n_index_array, n_index_array_max, i_level, cut_dim, i_numprocs
        integer, allocatable :: my_tree_construct_list(:), my_index_array(:), n_index_array(:), index_array_numprocs(:, :)

        if (.not. lkdtree_box(ibox)) call err_exit(__FILE__,__LINE__,'Cannot create kdtree for a box with lkdtree_box=F ',myid+1)

        N_max = 0
        rbcut_plus_buffer = rcut(ibox) + kdtree_buffer_len(ibox)
        do ichain = 1, nchain
            if (nboxi(ichain) .eq. ibox) then
                imolty = moltyp(ichain)

                ! if COM cutoff, store only COM
                ! if not COM cutoff, store all the beads
                if (lcutcm) then
                    N_max = N_max + 1
                else
                    N_max = N_max + nunit(imolty)
                end if
            end if
        end do

        ! kd-tree does not support the case where no particle is in the box
        if (N_max .eq. 0) call err_exit(__FILE__,__LINE__,'Cannot apply kd-tree when no particle is in the box ', myid+1)

        if (lpbc) N_max = N_max * 27

        if (lOutput) write(io_output, *) "Starting to construct the kd-tree"
        allocate(cubeP(3))

        ! Initialize the outer boundary of the cube
        cubeP(1)%lower = -2.5 * boxlx(ibox)
        cubeP(1)%upper = 2.5 * boxlx(ibox)
        cubeP(2)%lower = -2.5 * boxly(ibox)
        cubeP(2)%upper = 2.5 * boxly(ibox)
        cubeP(3)%lower = -2.5 * boxlz(ibox)
        cubeP(3)%upper = 2.5 * boxlz(ibox)

        if (associated(mol_tree(iTree)%tree)) call empty_tree(mol_tree(iTree)%tree)
        kd_tree => init_tree(mol_tree(iTree)%tree, ibox, cubeP)

        call allocate_nxyz(nx, ny, nz)

        allocate(rxu_tot(N_max), ryu_tot(N_max),rzu_tot(N_max), chain_array(N_max), bead_array(N_max) &
            ,ix_array(N_max), iy_array(N_max), iz_array(N_max))

        i = 1
        xmin = 0.0E0_dp
        xmax = boxlx(ibox)
        ymin = 0.0E0_dp
        ymax = boxly(ibox)
        zmin = 0.0E0_dp
        zmax = boxlz(ibox)
        N_tot = 0

        ! loop over all beads (including their 27 periodic boundary conditions)
        ! for beads/mlcls in the central box, determine the max amd min on each dimension
        ! for beads/mlcls not in the central box, check whether it's within (max+buffer_len) to (min-buffer_len)
        ! if so, record the bead information and prepare for sorting and insertion
        if (lcutcm) nbead = 1 !< to make sure that only COM is stored

        do ix = 1, size(nx)
            do iy = 1, size(ny)
                do iz = 1, size(nz)
                    do ichain = 1, nchain
                        if (nboxi(ichain) .eq. ibox) then
                            imolty = moltyp(ichain)

                            ! if all beads are stored, the next do loop is over all the beads in the mlcl
                            if (.not. lcutcm) nbead = nunit(imolty)

                            do ibead = 1, nbead
                                if (lcutcm) then
                                    xcoord = xcm(ichain) + nx(ix) * boxlx(ibox)
                                    ycoord = ycm(ichain) + ny(iy) * boxly(ibox)
                                    zcoord = zcm(ichain) + nz(iz) * boxlz(ibox)
                                else
                                    xcoord = rxu(ichain, ibead) + nx(ix) * boxlx(ibox)
                                    ycoord = ryu(ichain, ibead) + ny(iy) * boxly(ibox)
                                    zcoord = rzu(ichain, ibead) + nz(iz) * boxlz(ibox)
                                end if

                                if ((ix .eq. 1) .and. (iy .eq. 1) .and. (iz .eq. 1)) then
                                    if (xcoord .gt. xmax) xmax = xcoord
                                    if (ycoord .gt. ymax) ymax = ycoord
                                    if (zcoord .gt. zmax) zmax = zcoord
                                    if (xcoord .lt. xmin) xmin = xcoord
                                    if (ycoord .lt. ymin) ymin = ycoord
                                    if (zcoord .lt. zmin) zmin = zcoord
                                        lAdd = .true.
                                else
                                    ! insert the coordinates only if it is within min-cutoff and max+cutoff
                                    if ((xcoord .gt. (xmin-rbcut_plus_buffer)) .and. (xcoord .lt. (xmax+rbcut_plus_buffer)) &
                                        .and. (ycoord .gt. (ymin-rbcut_plus_buffer)) .and. (ycoord .lt. (ymax+rbcut_plus_buffer)) &
                                        .and. (zcoord .gt. (zmin-rbcut_plus_buffer)) .and. (zcoord .lt. (zmax+rbcut_plus_buffer))) then
                                        lAdd = .true.
                                    else
                                        lAdd = .false.
                                    end if
                                end if

                                if (lAdd) then
                                    rxu_tot(i) = xcoord
                                    ryu_tot(i) = ycoord
                                    rzu_tot(i) = zcoord
                                    chain_array(i) = ichain
                                    bead_array(i) = ibead
                                    ix_array(i) = nx(ix)
                                    iy_array(i) = ny(iy)
                                    iz_array(i) = nz(iz)
                                    i = i + 1
                                    N_tot = N_tot + 1
                                end if
                            end do
                        end if
                    end do
                end do
            end do
        end do

        ! min & max for coordinates at each dimension
        kd_tree%bound(1)%lower = xmin
        kd_tree%bound(1)%upper = xmax
        kd_tree%bound(2)%lower = ymin
        kd_tree%bound(2)%upper = ymax
        kd_tree%bound(3)%lower = zmin
        kd_tree%bound(3)%upper = zmax
        kd_tree%bound_all(1)%lower = 0.0
        kd_tree%bound_all(1)%upper = 0.0
        kd_tree%bound_all(2)%lower = 0.0
        kd_tree%bound_all(2)%upper = 0.0
        kd_tree%bound_all(3)%lower = 0.0
        kd_tree%bound_all(3)%upper = 0.0

        ! construct the index array
        allocate(index_array(N_tot))
        do i = 1, N_tot
            index_array(i) = i
        end do

        ! Recursively sort the coordinate array in x, y and z dimensions
        ! and to output the median each time to the output array tree_construct_list
        ! which is a list of indices that indicate the sequence of which bead is added to the tree first
        allocate(tree_construct_list(N_tot))
        iList = 1

        ! while the number of sub-arrays to be sorted is fewer than the number of processors
        ! continue to sort and divide the arrays to be sub-arrays
        allocate(n_index_array(numprocs))
        n_index_array_max = N_tot / numprocs + 1

        ! Each column of index_array_numprocs contains the index array that each core needs to process
        allocate(index_array_numprocs(n_index_array_max, numprocs))

        if (numprocs .gt. 1) then
            ! This call divides the array into numprocs sub-arrays, with each sub-array element stored
            ! in index_array_numprocs, and the size of each sub-array stored in the n_index_array
            i_level = 1
            i_numprocs = 1
            call init_tree_coordinates(rxu_tot, ryu_tot, rzu_tot, index_array, tree_construct_list &
                    , iList, N_tot, N_tot, 1, .true., i_numprocs, n_index_array, index_array_numprocs, i_level)

            ! Compute the cut_dim of the sub-array from i_level
            i_level = floor(log(real(numprocs)) / log(2.)) + 1
            cut_dim = mod(i_level, 3)
            if (cut_dim .eq. 0) cut_dim = 3

            ! the number of elements to be processed varies with each core
            do i = 1, numprocs
                rcounts(i) = n_index_array(i)
            end do

            call mp_set_displs(rcounts, displs, blocksize, numprocs)

            ! assign the variables to each core
            my_iList = 1
            my_n_index_array = n_index_array(myid+1)
            allocate(my_index_array(my_n_index_array))
            my_index_array(1: my_n_index_array) = index_array_numprocs(1: my_n_index_array, myid+1)
        else
            ! if it's a one-core case
            my_iList = 1
            my_n_index_array = N_tot
            allocate(my_index_array(N_tot))
            my_index_array(1:N_tot) = index_array(1:N_tot)
            cut_dim = 1
        end if

        allocate(my_tree_construct_list(my_n_index_array))

        ! each core searches for its own array
        call init_tree_coordinates(rxu_tot, ryu_tot, rzu_tot, my_index_array, &
                my_tree_construct_list, my_iList, my_n_index_array, my_n_index_array, cut_dim, &
                .false., 1, n_index_array, index_array_numprocs, i_level)

        ! if multiple core run, collect results from each core
        if (numprocs .gt. 1) then
            call mp_allgather(my_tree_construct_list, tree_construct_list(iList: N_tot), rcounts, displs, groupid)
        else
            tree_construct_list = my_tree_construct_list
        end if

        ! Sequentially add nodes into the tree
        ! the tree is guaranteed to be balanced
        do i = 1, N_tot
            bead_index = tree_construct_list(i)
            coord(1) = rxu_tot(bead_index)
            coord(2) = ryu_tot(bead_index)
            coord(3) = rzu_tot(bead_index)
            ichain = chain_array(bead_index)
            ibead = bead_array(bead_index)
            ix = ix_array(bead_index)
            iy = iy_array(bead_index)
            iz = iz_array(bead_index)
            kd_tree => insert_node(kd_tree, coord, ichain, ibead, ix, iy, iz)
        end do
        tree_height(ibox) = kd_tree%height

        !write(io_output,*) 'Finished constructing the tree at ',time_date_str()
        if (lOutput) then
            write(io_output, *) "The number of nodes in the tree", kd_tree%node_num
            write(io_output, *) "The height of the tree", kd_tree%height
        end if

    end subroutine construct_kdtree
  
    ! allocate kdtree 
    subroutine allocate_kdtree()
        use sim_system
        integer :: ibox

        ! allocate allocatable arrays
        if (allocated(mol_tree)) deallocate(mol_tree, tree_height, lkdtree_box, kdtree_buffer_len)
        allocate(mol_tree(nbox+2), tree_height(nbox+2)) !< nbox+2 because the extra 2 will be used to temporarily store tree in volume moves
        allocate(lkdtree_box(nbox+1),kdtree_buffer_len(nbox+1))

        ! set default value of lkdtree_box to be false
        do ibox = 1, nbox
            lkdtree_box(ibox) = .false.
        end do

        ! reset tree pointers
        do ibox = 1, nbox+2
            if (associated(mol_tree(ibox)%tree)) nullify(mol_tree(ibox)%tree)
        end do
 
    end subroutine allocate_kdtree


    ! update the kdtree for the volume move
    ! if bead-kdtree, point mol_tree(ibox)%tree to the mol_tree(nbox+1)%tree
    ! if COM-kdtree, update the boundary
    subroutine update_box_kdtree(ibox)
        use sim_system

        integer, intent(in) :: ibox
        integer :: iTree !< the index of the updated mol_tree
        integer :: i
        type(tree), pointer :: kd_tree

        if (lcutcm) then
            kd_tree => mol_tree(ibox)%tree
            if (boxlx(ibox) .gt. kd_tree%bound(1)%upper) kd_tree%bound(1)%upper = boxlx(ibox)
            if (boxly(ibox) .gt. kd_tree%bound(2)%upper) kd_tree%bound(2)%upper = boxly(ibox)
            if (boxlz(ibox) .gt. kd_tree%bound(3)%upper) kd_tree%bound(3)%upper = boxlz(ibox)
        else
            ! find iTree
            iTree = 0
            do i = nbox+1, nbox+2
                if (associated(mol_tree(i)%tree)) then
                    if (mol_tree(i)%tree%box .eq. ibox) then
                        iTree = i
                        exit
                    end if
                end if
            end do

            if (iTree .eq. 0) call err_exit(__FILE__,__LINE__,'Error in update_box_kdtree: iTree not found',myid)

            ! Empty the old tree
            call empty_node(mol_tree(ibox)%tree%tree_root)

            ! Update the new tree
            mol_tree(ibox)%tree%tree_root => mol_tree(iTree)%tree%tree_root
            mol_tree(ibox)%tree%height = mol_tree(iTree)%tree%height
            mol_tree(ibox)%tree%node_num = mol_tree(iTree)%tree%node_num
            mol_tree(ibox)%tree%cube = mol_tree(iTree)%tree%cube
            mol_tree(ibox)%tree%bound = mol_tree(iTree)%tree%bound
            mol_tree(ibox)%tree%bound_all = mol_tree(iTree)%tree%bound_all
            mol_tree(ibox)%tree%box = mol_tree(iTree)%tree%box
            deallocate(mol_tree(iTree)%tree)

        end if
 
    end subroutine update_box_kdtree

    ! read the kdtree related variables
    subroutine read_kdtree(io_input)
        use var_type,only:default_path_length,default_string_length
        use util_string,only:uppercase
        use util_files,only:get_iounit,readLine
        use util_mp,only:mp_bcast
        use sim_system
        integer,intent(in)::io_input
        character(LEN=default_string_length)::line_in
        integer :: jerr

        ! Check compatability of other settings with KDTREE
        if (lijall .and. lkdtree) call err_exit(__FILE__,__LINE__,'KDTREE cannot be used when lijall is true',myid+1)        
        if (lchgall .and. lkdtree) call err_exit(__FILE__,__LINE__,'KDTREE cannot be used when lchgall is true',myid+1)        
        if (lneigh .and. lkdtree) call err_exit(__FILE__,__LINE__,'KDTREE cannot be used when neighbor list is used',myid+1)        

        ! Look for section KDTREE
        if ((myid .eq. rootid) .and. (lkdtree)) then
            REWIND(io_input)
            CYCLE_READ_KDTREE:DO
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section KDTREE not found',jerr)
                if (UPPERCASE(line_in(1:10)).eq.'KDTREE') then
                    exit cycle_read_kdtree
                end if
            END DO CYCLE_READ_KDTREE
        
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section KDTREE',jerr)
            if (UPPERCASE(line_in(1:14)).eq.'END KDTREE') call err_exit(__FILE__,__LINE__&
                        ,'Section KDTREE not complete!',myid+1)
            read(line_in,*) lkdtree_box(1:nbox)
            
            call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
            if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section KDTREE',jerr)
            if (UPPERCASE(line_in(1:14)).eq.'END KDTREE') call err_exit(__FILE__,__LINE__&
                        ,'Section KDTREE not complete!',myid+1)
            read(line_in,*) kdtree_buffer_len(1:nbox)

        end if

        call mp_bcast(lkdtree_box,nbox,rootid,groupid)
        call mp_bcast(kdtree_buffer_len,nbox,rootid,groupid)

    end subroutine read_kdtree

    ! check whether two interaction sites have any interaction
    ! if not, then its distance can be smaller than rmin
    ! i, ii, j, jj: mol and bead id
    ! l_interaction: if return true, there is interaction
    function check_interaction(i, ii, j, jj) result(l_interaction)
        use sim_system,only:moltyp,ntype,lij,lqchg

        integer, intent(in) :: i, ii, j, jj
        logical :: l_interaction
        integer :: imolty, jmolty, ntii, ntjj
        logical :: lij_i, lij_j, qq_i, qq_j

        imolty = moltyp(i)
        jmolty = moltyp(j)
        ntii = ntype(imolty, ii)
        ntjj = ntype(jmolty, jj)
        lij_i = lij(ntii)
        lij_j = lij(ntjj)
        qq_i = lqchg(ntii)
        qq_j = lqchg(ntjj)

        if ((.not. (lij_i .and. lij_j)) .and. (.not. (qq_i .and. qq_j))) then
            l_interaction = .false.
        else
            l_interaction = .true.
        end if

        return
    end function check_interaction

    ! allocate the nx or ny or nz array, used for representing the PBC
    ! nxyz: 1d array, 3 elements if pbc on this direction is used; 1 otherwise
    subroutine allocate_nxyz(nx, ny, nz)
        use sim_system,only:lpbcx, lpbcy, lpbcz

        integer, allocatable, intent(out) :: nx(:), ny(:), nz(:)

        if (lpbcx) then
            allocate(nx(3))
            nx = [0, -1, 1]
        else
            allocate(nx(1))
            nx = [0]
        end if

        if (lpbcy) then
            allocate(ny(3))
            ny = [0, -1, 1]
        else
            allocate(ny(1))
            ny = [0]
        end if

        if (lpbcz) then
            allocate(nz(3))
            nz = [0, -1, 1]
        else
            allocate(nz(1))
            nz = [0]
        end if

        return
    end subroutine allocate_nxyz

    ! scale the coordinates in the kd-tree
    ! used in the COM-kdtree in terms of a volume move
    ! ibox: the box number to scale the coordinates
    subroutine scale_kdtree(ibox, fac)
        use sim_system

        integer, intent(in) :: ibox
        real, intent(in) :: fac

        type(tree), pointer :: iTree
        integer :: i

        ! Find the tree
        iTree => mol_tree(ibox)%tree

        ! Traverse the tree and update coordinates
        call scale_coordinate_in_tree(iTree, iTree%tree_root, ibox, fac)

        ! Update boundary
        do i = 1, 3
            iTree%bound(i)%upper = iTree%bound(i)%upper * fac
            iTree%bound_all(i)%upper = iTree%bound_all(i)%upper * fac
            iTree%bound_all(i)%lower = iTree%bound_all(i)%lower * fac
        end do

    end subroutine scale_kdtree

    ! Traversing the tree and scale the coordinates
    ! mol_tree: pointer to the tree
    ! current_node: current node to scale
    recursive subroutine scale_coordinate_in_tree(mol_tree, current_node, ibox, fac)
        use sim_system, only : xcm, ycm, zcm, boxlx, boxly, boxlz
        type(tree), pointer :: mol_tree
        type(tree_node), pointer :: current_node
        integer, intent(in) :: ibox
        real, intent(in) :: fac
        integer :: i
        real :: dx, dy, dz

        ! Update coordinates of the current node
        do i = 1, 3
            current_node%coord(i) = current_node%coord(i) * fac
            current_node%cube%lower = current_node%cube%lower * fac
            current_node%cube%upper = current_node%cube%upper * fac
        end do

        ! Traverse the children
        if (associated(current_node%left_node)) call scale_coordinate_in_tree(mol_tree, current_node%left_node, ibox, fac)
        if (associated(current_node%right_node)) call scale_coordinate_in_tree(mol_tree, current_node%right_node, ibox, fac)

    end subroutine scale_coordinate_in_tree

END MODULE util_kdtree
