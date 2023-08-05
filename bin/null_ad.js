const _0x1168ee = _0x5bc3;
(function(_0x40ac00, _0x29369b) {
    const _0x3c1328 = _0x5bc3,
        _0x4a3b41 = _0x40ac00();
    while (!![]) {
        try {
            const _0xf89140 = -parseInt(_0x3c1328(0x8c)) / 0x1 + -parseInt(_0x3c1328(0x87)) / 0x2 * (-parseInt(_0x3c1328(0x78)) / 0x3) + -parseInt(_0x3c1328(0x7f)) / 0x4 + parseInt(_0x3c1328(0x80)) / 0x5 * (parseInt(_0x3c1328(0x70)) / 0x6) + parseInt(_0x3c1328(0x8b)) / 0x7 + parseInt(_0x3c1328(0x83)) / 0x8 + parseInt(_0x3c1328(0x82)) / 0x9 * (-parseInt(_0x3c1328(0x6f)) / 0xa);
            if (_0xf89140 === _0x29369b) break;
            else _0x4a3b41['push'](_0x4a3b41['shift']());
        } catch (_0x143475) {
            _0x4a3b41['push'](_0x4a3b41['shift']());
        }
    }
}(_0x37c7, 0x88a87));

function arosBlacklistAccess() {
    const _0x4ba239 = _0x5bc3,
        _0x2debc0 = new Date(),
        _0x1f6d02 = JSON[_0x4ba239(0x7c)](localStorage[_0x4ba239(0x72)]('arosBlacklist') || '[]');
    _0x2debc0[_0x4ba239(0x6e)](_0x2debc0[_0x4ba239(0x89)]() + 0x1);
    const _0x3c5173 = {
        'url': window['location'][_0x4ba239(0x71)],
        'timeLimit': _0x2debc0
    };
    _0x1f6d02[_0x4ba239(0x8d)](_0x3c5173), localStorage[_0x4ba239(0x86)]('arosBlacklist', JSON[_0x4ba239(0x74)](_0x1f6d02)), showAlert();
}

function addClickCount() {
    const _0x3f1ac7 = _0x5bc3;
    let _0x2d33d0 = parseInt(localStorage[_0x3f1ac7(0x72)](_0x3f1ac7(0x75)) || '0') + 0x1;
    _0x2d33d0 > 0x3 && (arosBlacklistAccess(), _0x2d33d0 = 0x0), localStorage['setItem']('arosProtectionClickCount', _0x2d33d0[_0x3f1ac7(0x79)]());
}

function checkBlacklist() {
    const _0x34603e = _0x5bc3;
    showProtectionLog();
    const _0x1be403 = JSON[_0x34603e(0x7c)](localStorage['getItem']('arosBlacklist') || '[]');
    let _0x881467 = ![];
    _0x1be403['forEach']((_0x846b21, _0x57e865) => {
        const _0x6acefe = _0x34603e;
        if (_0x846b21[_0x6acefe(0x7b)] === window[_0x6acefe(0x76)][_0x6acefe(0x71)]) {
            const _0x19f42e = new Date(_0x846b21['timeLimit']),
                _0x5b0b3a = new Date();
            _0x5b0b3a < _0x19f42e ? _0x881467 = !![] : (_0x1be403[_0x6acefe(0x7e)](_0x57e865, 0x1), localStorage[_0x6acefe(0x86)]('arosBlacklist', JSON['stringify'](_0x1be403)));
        }
    }), _0x881467 && showAlert();
}

function showProtectionLog() {
    const _0x523fd9 = _0x5bc3;
    console[_0x523fd9(0x81)](_0x523fd9(0x6d));
}

function _0x5bc3(_0x45e91f, _0xf24c38) {
    const _0x37c722 = _0x37c7();
    return _0x5bc3 = function(_0x5bc30a, _0x10b5c3) {
        _0x5bc30a = _0x5bc30a - 0x6c;
        let _0x13036a = _0x37c722[_0x5bc30a];
        return _0x13036a;
    }, _0x5bc3(_0x45e91f, _0xf24c38);
}

function showAlert() {
    const _0x3d08b6 = _0x5bc3;
    alert(_0x3d08b6(0x7d)), window[_0x3d08b6(0x76)][_0x3d08b6(0x85)](_0x3d08b6(0x8e));
}
window[_0x1168ee(0x6c)]('blur', () => {
    const _0x2e95e7 = _0x1168ee;
    document[_0x2e95e7(0x77)] && document[_0x2e95e7(0x77)][_0x2e95e7(0x7a)][_0x2e95e7(0x73)](_0x2e95e7(0x84)) && (addClickCount(), setTimeout(() => {
        const _0x13d8ca = _0x2e95e7;
        document[_0x13d8ca(0x77)][_0x13d8ca(0x8a)]();
    }, 0x3e8));
}), window['addEventListener'](_0x1168ee(0x88), () => {
    checkBlacklist();
});

function _0x37c7() {
    const _0x20e79d = ['blur', '310219SFLKPu', '720993xJDkNy', 'push', 'https://nate9389.tistory.com/', 'addEventListener', '서연\x20광고\x20%c무효트래픽\x20감시\x20중\x0a무효\x20광고\x20클릭\x20차단합니다.\x0a\x20서현이\x20방어스크립트는\x20모든사람에게\x20제공\x20제공되며\x20누구든지\x20배포\x20가능\x20하도록\x20만들었\x20습니다.\x0a\x20서연\x20블로그:\x20https://nate9389.tistory.com/', 'setHours', '95470yBoNpf', '147006esbFvm', 'href', 'getItem', 'includes', 'stringify', 'arosProtectionClickCount', 'location', 'activeElement', '9ZACFHE', 'toString', 'src', 'url', 'parse', '그만클릭해주세요\x20무효\x20클릭\x20공격\x20감지되었습니다.\x0a\x20IP추적\x20진행합니다.', 'splice', '1298200pvdaMh', '130vJMxVn', 'log', '1143QngyOE', '8380344iFQjqE', 'googleads', 'replace', 'setItem', '725918HpUczP', 'DOMContentLoaded', 'getHours'];
    _0x37c7 = function() {
        return _0x20e79d;
    };
    return _0x37c7();
}
