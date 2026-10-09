#include "Electrospray/IO/ElectrosprayRunSession.h"

#include <cstdio>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <io.h>
#else
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#endif

namespace
{
namespace fs = std::filesystem;

void rejectRedirectedPath(const fs::path& path)
{
    // Windows junctions are reparse points too, not just directory symlinks.
#ifdef _WIN32
    const DWORD attributes = GetFileAttributesW(path.c_str());
    if (attributes != INVALID_FILE_ATTRIBUTES && (attributes & FILE_ATTRIBUTE_REPARSE_POINT))
        throw std::runtime_error("拒绝清理链接或重解析路径：" + path.string());
#endif
    if (fs::is_symlink(fs::symlink_status(path)))
        throw std::runtime_error("拒绝清理符号链接：" + path.string());
}

struct OutputInventory
{
    bool hasEntries = false;
    std::uintmax_t files = 0;
};

OutputInventory inspectOutput(const fs::path& root)
{
    const auto path = root / "output";
    rejectRedirectedPath(path);
    if (!fs::exists(path)) return {};
    if (!fs::is_directory(path) || fs::canonical(path).parent_path() != root)
        throw std::runtime_error("输出路径必须是当前算例内的普通 output 目录：" + path.string());
    OutputInventory result;
    // Validate the WHOLE tree before removing anything; never follow links.
    for (const auto& entry : fs::recursive_directory_iterator(path))
    {
        result.hasEntries = true;
        rejectRedirectedPath(entry.path());
        if (entry.is_regular_file()) ++result.files;
        else if (!entry.is_directory())
            throw std::runtime_error("输出目录含非普通文件，未清理旧结果：" + entry.path().string());
    }
    return result;
}
}

struct ElectrosprayRunSession::Impl
{
    std::filesystem::path root;
#ifdef _WIN32
    HANDLE lock = INVALID_HANDLE_VALUE;
#else
    int lock = -1;
#endif

    ~Impl()
    {
#ifdef _WIN32
        if (lock != INVALID_HANDLE_VALUE) CloseHandle(lock); // delete-on-close, including process termination
#else
        // Keep the inode: unlinking a flock file could let a third process lock
        // a different inode while a second process already holds the old one.
        if (lock >= 0) { flock(lock, LOCK_UN); close(lock); }
#endif
    }
};

ElectrosprayRunSession::ElectrosprayRunSession(const fs::path& caseDirectory)
    : impl_(std::make_unique<Impl>())
{
    impl_->root = fs::canonical(caseDirectory);
    if (impl_->root == impl_->root.root_path()
        || !fs::is_regular_file(impl_->root / "system" / "run.dict"))
        throw std::runtime_error("启动目录必须是包含 system/run.dict 的算例目录");
    const auto path = impl_->root / ".electrospray.lock";
    rejectRedirectedPath(path);
    if (fs::exists(path) && (!fs::is_regular_file(path) || fs::hard_link_count(path) != 1))
        throw std::runtime_error("算例运行锁路径不是普通独立文件：" + path.string());
#ifdef _WIN32
    impl_->lock = CreateFileW(path.c_str(), GENERIC_READ | GENERIC_WRITE | DELETE, 0, nullptr,
        OPEN_ALWAYS, FILE_ATTRIBUTE_HIDDEN | FILE_FLAG_DELETE_ON_CLOSE | FILE_FLAG_OPEN_REPARSE_POINT, nullptr);
    if (impl_->lock == INVALID_HANDLE_VALUE)
        throw std::runtime_error("无法独占当前算例；可能已有计算正在运行，或目录不可写。未清理旧结果：" + impl_->root.string());
#else
    impl_->lock = open(path.c_str(), O_CREAT | O_RDWR | O_NOFOLLOW, 0600);
    if (impl_->lock < 0 || flock(impl_->lock, LOCK_EX | LOCK_NB) != 0)
        throw std::runtime_error("无法独占当前算例；可能已有计算正在运行，或目录不可写。未清理旧结果：" + impl_->root.string());
#endif
}

ElectrosprayRunSession::~ElectrosprayRunSession() = default;

bool ElectrosprayRunSession::terminalInput()
{
#ifdef _WIN32
    return _isatty(_fileno(stdin)) != 0;
#else
    return isatty(fileno(stdin)) != 0;
#endif
}

bool ElectrosprayRunSession::prepare(RunStartMode mode, bool interactive,
    std::istream& input, std::ostream& output)
{
    const auto inventory = inspectOutput(impl_->root);
    if (mode == RunStartMode::Automatic)
    {
        if (interactive) mode = RunStartMode::Ask;
        else if (inventory.hasEntries)
            throw std::runtime_error("当前算例已有输出，非交互运行未获覆盖确认；请使用 --start new 从初态重算并覆盖 output/。旧结果未改动。");
        else mode = RunStartMode::New;
    }
    if (mode == RunStartMode::Ask)
    {
        output << "[电喷雾] 请选择启动方式：\n"
            << "  1  从初态重新计算（覆盖当前算例 output/ 内全部旧结果，不可由程序恢复）\n"
            << "  2  接上次计算（断点续算尚未支持）\n"
            << "  0  取消，保留旧结果\n"
            << "  当前输出目录：" << (impl_->root / "output").string()
            << "，已有文件数=" << inventory.files << '\n';
        for (;;)
        {
            output << "请选择 [0/1/2，回车取消]：" << std::flush;
            std::string answer;
            if (!std::getline(input, answer) || answer.find_first_not_of(" \t\r") == std::string::npos)
                answer = "0";
            else
            {
                answer.erase(0, answer.find_first_not_of(" \t\r"));
                answer.erase(answer.find_last_not_of(" \t\r") + 1);
            }
            if (answer == "0")
            {
                output << "[电喷雾] 已取消启动；旧结果未改动。\n";
                return false;
            }
            if (answer == "1") { mode = RunStartMode::New; break; }
            if (answer == "2") { mode = RunStartMode::Resume; break; }
            output << "请输入 0、1 或 2。\n";
        }
    }
    if (mode == RunStartMode::Resume)
        throw std::runtime_error("电喷雾断点续算尚未支持；现有 DAT/CSV 不含完整恢复状态。本次未启动，旧结果未改动。");
    if (mode != RunStartMode::New) throw std::invalid_argument("未知的电喷雾启动方式");

    const auto target = impl_->root / "output";
    const auto confirmed = inspectOutput(impl_->root); // revalidate after a possibly long prompt
    if (fs::exists(target)) fs::remove_all(target); // exact validated case-local directory, never cache/system
    output << "[电喷雾] 从初态重新计算；已清理旧输出文件 " << confirmed.files
        << " 个，目录：" << target.string() << "。system/ 和 cache/ 保留。\n" << std::flush;
    return true;
}
