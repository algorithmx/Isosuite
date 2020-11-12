global const default_settings_findsym = [
    "monoclinicAxes" => "a(b)c",
    #Axes of a monoclinic space group.  Enter a(b)c, c(-b)a, ab(c), ba(-c), (a)bc, 
    #or (-a)cb.  The axis in parentheses is the unique axis.  
    #The default value is a(b)c.
    "monoclinicCell" => "1",
    # Cell choice of a monoclinic space group.  Enter 1, 2, or 3.  The default value is 1.
    "orthorhombicAxes" => "abc",
    #Axes of an orthorhombic space group.  Enter abc, ba-c, cab, -cba, bca, or a-cb.
    #The default value is abc.
    "originChoice" => "2",
    #Origin choice for a space group.  Enter 1 or 2.  The default value is 2.
    "hexagonalAxes" => "",
    #Use hexagonal axes for R-centered hexagonal space groups.  This is the default value.
    #"rhombohedralAxes" => "",
    #Use rhombohedral axes for R-centered hexagonal space groups.  The default is to
    #use hexagonal axes.
]

# -------------------------------------------------
#* QE PW reference
#* https://www.quantum-espresso.org/Doc/INPUT_PW.html#rhombohedral

global const QE_default_equivalent_settings_findsym = [
    "monoclinicAxes" => "ab(c)", 
    #* uniqueb default is FALSE, the two fold axis or the mirror normal is parallel to the c axis.
    #Axes of a monoclinic space group.  Enter a(b)c, c(-b)a, ab(c), ba(-c), (a)bc, 
    #or (-a)cb.  The axis in parentheses is the unique axis. 
    #The default value is a(b)c.
    "monoclinicCell" => "1",
    # Cell choice of a monoclinic space group.  Enter 1, 2, or 3.  The default value is 1.
    "orthorhombicAxes" => "abc",
    #Axes of an orthorhombic space group.  Enter abc, ba-c, cab, -cba, bca, or a-cb.
    #The default value is abc.
    "originChoice" => "1",
    #* origin_choice default is 1
    #Origin choice for a space group.  Enter 1 or 2.  The default value is 2.
    #"hexagonalAxes" => "",
    #Use hexagonal axes for R-centered hexagonal space groups.  This is the default value.
    "rhombohedralAxes" => "",
    #* rhombohedral default is TRUE
    #Use rhombohedral axes for R-centered hexagonal space groups.  The default is to
    #use hexagonal axes.
]

