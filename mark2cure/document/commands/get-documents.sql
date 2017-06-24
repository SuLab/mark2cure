SELECT  `document_document`.`id`,
        `document_document`.`document_id` as `pmid`,
        `document_section`.`kind` as `section`,
        `document_section`.`id` as `section_id`,
        `document_section`.`text`

FROM `document_document`

LEFT JOIN `document_section`
    ON `document_section`.`document_id` = `document_document`.`id`

WHERE `document_document`.`id` IN ({0})

