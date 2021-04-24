import numpy as np
cimport numpy as np

cdef extern from 'wrapped_code_23.h':
    void autofunc(double _Dummy_112, double _Dummy_113, double _Dummy_114, double _Dummy_115, double _Dummy_116, double _Dummy_117, double _Dummy_118, double _Dummy_119, double _Dummy_120, double _Dummy_121, double _Dummy_122, double _Dummy_123, double _Dummy_124, double _Dummy_125, double _Dummy_126, double _Dummy_22, double _Dummy_23, double _Dummy_24, double _Dummy_25, double _Dummy_26, double _Dummy_27, double _Dummy_28, double _Dummy_29, double _Dummy_30, double _Dummy_31, double _Dummy_32, double _Dummy_33, double _Dummy_34, double _Dummy_35, double _Dummy_36, double _Dummy_37, double _Dummy_38, double _Dummy_39, double _Dummy_40, double _Dummy_41, double _Dummy_42, double _Dummy_43, double _Dummy_44, double _Dummy_45, double _Dummy_46, double _Dummy_47, double _Dummy_48, double _Dummy_49, double _Dummy_50, double _Dummy_51, double _Dummy_52, double _Dummy_53, double _Dummy_54, double _Dummy_55, double _Dummy_56, double _Dummy_57, double _Dummy_58, double _Dummy_59, double _Dummy_60, double _Dummy_61, double _Dummy_62, double _Dummy_63, double _Dummy_64, double _Dummy_65, double _Dummy_66, double _Dummy_67, double _Dummy_68, double _Dummy_69, double _Dummy_70, double _Dummy_71, double _Dummy_72, double _Dummy_73, double _Dummy_74, double _Dummy_75, double _Dummy_76, double _Dummy_77, double _Dummy_78, double _Dummy_79, double _Dummy_80, double _Dummy_81, double _Dummy_82, double _Dummy_83, double _Dummy_84, double _Dummy_85, double _Dummy_86, double _Dummy_87, double _Dummy_88, double _Dummy_89, double _Dummy_90, double _Dummy_91, double _Dummy_92, double _Dummy_93, double _Dummy_94, double _Dummy_95, double _Dummy_96, double _Dummy_97, double _Dummy_98, double _Dummy_99, double _Dummy_100, double _Dummy_101, double _Dummy_102, double _Dummy_103, double _Dummy_104, double _Dummy_105, double _Dummy_106, double _Dummy_107, double _Dummy_108, double _Dummy_109, double _Dummy_110, double _Dummy_111, double *out_3473116775766001580)

def autofunc_c(double _Dummy_112, double _Dummy_113, double _Dummy_114, double _Dummy_115, double _Dummy_116, double _Dummy_117, double _Dummy_118, double _Dummy_119, double _Dummy_120, double _Dummy_121, double _Dummy_122, double _Dummy_123, double _Dummy_124, double _Dummy_125, double _Dummy_126, double _Dummy_22, double _Dummy_23, double _Dummy_24, double _Dummy_25, double _Dummy_26, double _Dummy_27, double _Dummy_28, double _Dummy_29, double _Dummy_30, double _Dummy_31, double _Dummy_32, double _Dummy_33, double _Dummy_34, double _Dummy_35, double _Dummy_36, double _Dummy_37, double _Dummy_38, double _Dummy_39, double _Dummy_40, double _Dummy_41, double _Dummy_42, double _Dummy_43, double _Dummy_44, double _Dummy_45, double _Dummy_46, double _Dummy_47, double _Dummy_48, double _Dummy_49, double _Dummy_50, double _Dummy_51, double _Dummy_52, double _Dummy_53, double _Dummy_54, double _Dummy_55, double _Dummy_56, double _Dummy_57, double _Dummy_58, double _Dummy_59, double _Dummy_60, double _Dummy_61, double _Dummy_62, double _Dummy_63, double _Dummy_64, double _Dummy_65, double _Dummy_66, double _Dummy_67, double _Dummy_68, double _Dummy_69, double _Dummy_70, double _Dummy_71, double _Dummy_72, double _Dummy_73, double _Dummy_74, double _Dummy_75, double _Dummy_76, double _Dummy_77, double _Dummy_78, double _Dummy_79, double _Dummy_80, double _Dummy_81, double _Dummy_82, double _Dummy_83, double _Dummy_84, double _Dummy_85, double _Dummy_86, double _Dummy_87, double _Dummy_88, double _Dummy_89, double _Dummy_90, double _Dummy_91, double _Dummy_92, double _Dummy_93, double _Dummy_94, double _Dummy_95, double _Dummy_96, double _Dummy_97, double _Dummy_98, double _Dummy_99, double _Dummy_100, double _Dummy_101, double _Dummy_102, double _Dummy_103, double _Dummy_104, double _Dummy_105, double _Dummy_106, double _Dummy_107, double _Dummy_108, double _Dummy_109, double _Dummy_110, double _Dummy_111):

    cdef np.ndarray[np.double_t, ndim=2] out_3473116775766001580 = np.empty((3,3))
    autofunc(_Dummy_112, _Dummy_113, _Dummy_114, _Dummy_115, _Dummy_116, _Dummy_117, _Dummy_118, _Dummy_119, _Dummy_120, _Dummy_121, _Dummy_122, _Dummy_123, _Dummy_124, _Dummy_125, _Dummy_126, _Dummy_22, _Dummy_23, _Dummy_24, _Dummy_25, _Dummy_26, _Dummy_27, _Dummy_28, _Dummy_29, _Dummy_30, _Dummy_31, _Dummy_32, _Dummy_33, _Dummy_34, _Dummy_35, _Dummy_36, _Dummy_37, _Dummy_38, _Dummy_39, _Dummy_40, _Dummy_41, _Dummy_42, _Dummy_43, _Dummy_44, _Dummy_45, _Dummy_46, _Dummy_47, _Dummy_48, _Dummy_49, _Dummy_50, _Dummy_51, _Dummy_52, _Dummy_53, _Dummy_54, _Dummy_55, _Dummy_56, _Dummy_57, _Dummy_58, _Dummy_59, _Dummy_60, _Dummy_61, _Dummy_62, _Dummy_63, _Dummy_64, _Dummy_65, _Dummy_66, _Dummy_67, _Dummy_68, _Dummy_69, _Dummy_70, _Dummy_71, _Dummy_72, _Dummy_73, _Dummy_74, _Dummy_75, _Dummy_76, _Dummy_77, _Dummy_78, _Dummy_79, _Dummy_80, _Dummy_81, _Dummy_82, _Dummy_83, _Dummy_84, _Dummy_85, _Dummy_86, _Dummy_87, _Dummy_88, _Dummy_89, _Dummy_90, _Dummy_91, _Dummy_92, _Dummy_93, _Dummy_94, _Dummy_95, _Dummy_96, _Dummy_97, _Dummy_98, _Dummy_99, _Dummy_100, _Dummy_101, _Dummy_102, _Dummy_103, _Dummy_104, _Dummy_105, _Dummy_106, _Dummy_107, _Dummy_108, _Dummy_109, _Dummy_110, _Dummy_111, <double*> out_3473116775766001580.data)
    return out_3473116775766001580