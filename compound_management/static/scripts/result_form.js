$(document).ready(function() {

     $('.db-radio').on('change', function() {
            console.log($(this).val()); // 선택된 체크박스의 값을 콘솔에 출력
     });

    $('input[name="dbOption"]').on('change', function() {
    const tableName = $(this).val();

    $.getJSON(`/get_columns/${tableName}/`, function(data) {
        const $table = $('#data-table');
        $table.empty();

//        const filteredFields = data.fields.filter(field => field.help_text.includes('check'));
//          const filteredFields = data.fields;
        const filteredFields = data.fields.filter(field => !('help_text' in field) || !field.help_text);
        console.log(filteredFields[1]);

        const $headerRow = $('<tr></tr>');
        filteredFields.forEach(function(field) {
            $headerRow.append(`<th>${field.name}</th>`);
        });
        $table.append($headerRow);

        const $inputRow = $('<tr></tr>');
        filteredFields.forEach(function(field) {
            $inputRow.append(`<td><input type="text" name="${field.name}"></td>`);
        });
        $table.append($inputRow);
    }
    );
});

});

function updateResultTable() {
    const $resultTableBody = $('#resultTable tbody');
    $resultTableBody.empty();

     const filteredFields = data.fields.filter(field => field.help_text.includes('check'));

    $('.data-checkbox:checked').each(function() {
        var $row = $(this).closest("tr")
        const data = $row.find("td[data-id]").data("id");
        const $newRow = $('<tr></tr>');
        const $cell = $('<td></td>').text(data);
        $newRow.append($cell);
        $resultTableBody.append($newRow);
    });
}


    // 페이지가 로드될 때, 저장된 탭 상태를 복원
    var activeTab = localStorage.getItem('activeTab');
    if (activeTab) {
        var $tabElement = $(activeTab);
        if ($tabElement.length) {
            $tabElement.click();
        }
    }

    // 탭 클릭 시, 로컬 스토리지에 상태 저장
    $('.tab-link').on('click', function() {
        var target = $(this).attr('href');
        localStorage.setItem('activeTab', target);
    });
$(document).on('change', '.data-checkbox', updateResultTable);