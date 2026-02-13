import Foundation
import BWACore

/// Owns an mmap'd file region. Read-only, MAP_PRIVATE.
/// Keeps the file descriptor and mapping alive via ARC.
public final class MappedFile: @unchecked Sendable {
    let pointer: UnsafeRawPointer
    let size: Int
    private let fd: Int32

    init(path: String) throws {
        fd = open(path, O_RDONLY)
        guard fd >= 0 else {
            throw BWAError.mmapFailed("Cannot open \(path)")
        }

        var st = stat()
        guard fstat(fd, &st) == 0 else {
            close(fd)
            throw BWAError.mmapFailed("Cannot stat \(path)")
        }
        size = Int(st.st_size)

        let ptr = mmap(nil, size, PROT_READ, MAP_PRIVATE, fd, 0)
        guard ptr != MAP_FAILED else {
            close(fd)
            throw BWAError.mmapFailed("mmap failed for \(path)")
        }
        pointer = UnsafeRawPointer(ptr!)
    }

    deinit {
        munmap(UnsafeMutableRawPointer(mutating: pointer), size)
        close(fd)
    }
}
