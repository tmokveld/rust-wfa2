pub(crate) struct CigarView<'a> {
    pub(crate) score: i32,
    pub(crate) end_v: i32,
    pub(crate) end_h: i32,
    pub(crate) operations: &'a [std::os::raw::c_char],
}

fn active_cigar_operations(
    operations: &[std::os::raw::c_char],
    begin_offset: i32,
    end_offset: i32,
) -> &[std::os::raw::c_char] {
    let Ok(begin_offset) = usize::try_from(begin_offset) else {
        return &[];
    };
    let Ok(end_offset) = usize::try_from(end_offset) else {
        return &[];
    };

    operations.get(begin_offset..end_offset).unwrap_or(&[])
}

fn operation_bytes(operations: &[std::os::raw::c_char]) -> &[u8] {
    // SAFETY: c_char is one byte, and the returned slice does not outlive `operations`.
    unsafe { std::slice::from_raw_parts(operations.as_ptr() as *const u8, operations.len()) }
}

pub(crate) fn swap_indel_ops_in_packed_cigar(cigar: &mut [u32]) {
    const SAM_CIGAR_OP_MASK: u32 = 0xF;
    const SAM_CIGAR_INS: u32 = 1;
    const SAM_CIGAR_DEL: u32 = 2;

    for encoded_op in cigar {
        let op_code = *encoded_op & SAM_CIGAR_OP_MASK;
        let swapped_op = match op_code {
            SAM_CIGAR_INS => SAM_CIGAR_DEL,
            SAM_CIGAR_DEL => SAM_CIGAR_INS,
            _ => continue,
        };
        *encoded_op = (*encoded_op & !SAM_CIGAR_OP_MASK) | swapped_op;
    }
}

pub(crate) fn swap_indel_ops_in_cigar_bytes(cigar: &mut [u8]) {
    for op in cigar {
        *op = match *op {
            b'I' => b'D',
            b'D' => b'I',
            other => other,
        };
    }
}

impl<'a> CigarView<'a> {
    pub(crate) fn new(
        score: i32,
        begin_offset: i32,
        end_offset: i32,
        end_v: i32,
        end_h: i32,
        operations: &'a [std::os::raw::c_char],
    ) -> Self {
        Self {
            score,
            end_v,
            end_h,
            operations: active_cigar_operations(operations, begin_offset, end_offset),
        }
    }

    pub(crate) fn active_operation_bytes(&self) -> &'a [u8] {
        operation_bytes(self.operations)
    }

    pub(crate) fn clipped_operations(&self, flank_len: usize) -> &[std::os::raw::c_char] {
        let Some(end_offset) = self.operations.len().checked_sub(flank_len) else {
            return &[];
        };

        if flank_len >= end_offset {
            return &[];
        }

        &self.operations[flank_len..end_offset]
    }

    #[cfg(test)]
    pub(crate) fn clipped_operation_bytes(&self, flank_len: usize) -> &[u8] {
        operation_bytes(self.clipped_operations(flank_len))
    }

    pub(crate) fn end_position(&self) -> Option<(usize, usize)> {
        if self.end_v < 0 || self.end_h < 0 {
            return None;
        }

        Some((self.end_v as usize, self.end_h as usize))
    }
}

pub type CigarOp = (usize, char);
